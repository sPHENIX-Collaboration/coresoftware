#include "ActsTransformations.h"
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <chrono>
using namespace std::chrono;

Acts::BoundSymMatrix ActsTransformations::rotateSvtxTrackCovToActs(
			        const SvtxTrack *track,
				Acts::GeometryContext geoCtxt) const
{
  Acts::BoundSymMatrix svtxCovariance = Acts::BoundSymMatrix::Zero();

  for(int i = 0; i < 6; ++i)
    {
      for(int j = 0; j < 6; ++j)
	{
	  svtxCovariance(i,j) = track->get_error(i,j);
	  /// Convert Svtx to mm and GeV units as Acts expects
	  if(i < 3 && j < 3)
	    svtxCovariance(i,j) *= Acts::UnitConstants::cm2;
	  else if (i < 3)
	    svtxCovariance(i,j) *= Acts::UnitConstants::cm;
	  else if (j < 3)
	    svtxCovariance(i,j) *= Acts::UnitConstants::cm;
	  
	}
    }

  printMatrix("svtx covariance, acts units: ", svtxCovariance);
 
  const double px = track->get_px();
  const double py = track->get_py();
  const double pz = track->get_pz();
  const double p = sqrt(px*px + py*py + pz*pz);
  
  const double uPx = px / p;
  const double uPy = py / p;
  const double uPz = pz / p;

  //Acts version
  const double cosTheta = uPz;
  const double sinTheta = sqrt(uPx * uPx + uPy * uPy);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = uPx * invSinTheta; // equivalent to x/r
  const double sinPhi = uPy * invSinTheta; // equivalent to y/r
    
  /// First we rotate to (x,y,z,time,Tx,Ty,Tz,q/p) to take advantage of the
  /// already created Acts rotation matrix from this basis into the Acts local basis
  /// We basically go backwards from rotateActsCovToSvtxTrack to get the Acts cov 
  /// from the SvtxTrack cov

  /// This is going from Acts->Svtx, so we will take the transpose
  Acts::ActsMatrixD<6,8> sphenixRot;
  sphenixRot.setZero();
  /// Make the xyz transform unity
  sphenixRot(0,0) = 1;
  sphenixRot(1,1) = 1;
  sphenixRot(2,2) = 1;
  sphenixRot(3,4) = 1./p;
  sphenixRot(4,5) = 1./p;
  sphenixRot(5,6) = 1./p;
  sphenixRot(3,7) = uPx * p * p;
  sphenixRot(4,7) = uPy * p * p;
  sphenixRot(5,7) = uPz * p * p;

  auto rotatedMatrix 
    = sphenixRot.transpose() * svtxCovariance * sphenixRot;
  
  /// Now take the 8x8 matrix and rotate it to Acts basis
  Acts::BoundToFreeMatrix jacobianLocalToGlobal = Acts::BoundToFreeMatrix::Zero();
  jacobianLocalToGlobal(0, Acts::eBoundLoc0) = -sinPhi;
  jacobianLocalToGlobal(0, Acts::eBoundLoc1) = -cosPhi * cosTheta;
  jacobianLocalToGlobal(1, Acts::eBoundLoc0) = cosPhi;
  jacobianLocalToGlobal(1, Acts::eBoundLoc1) = -sinPhi * cosTheta;
  jacobianLocalToGlobal(2, Acts::eBoundLoc1) = sinTheta;
  jacobianLocalToGlobal(3, Acts::eBoundTime) = 1;
  jacobianLocalToGlobal(4, Acts::eBoundPhi) = -sinTheta * sinPhi;
  jacobianLocalToGlobal(4, Acts::eBoundTheta) = cosTheta * cosPhi;
  jacobianLocalToGlobal(5, Acts::eBoundPhi) = sinTheta * cosPhi;
  jacobianLocalToGlobal(5, Acts::eBoundTheta) = cosTheta * sinPhi;
  jacobianLocalToGlobal(6, Acts::eBoundTheta) = -sinTheta;
  jacobianLocalToGlobal(7, Acts::eBoundQOverP) = 1;

  /// Since we are using the local to global jacobian we do R^TCR instead of 
  /// RCR^T
  Acts::BoundSymMatrix actsLocalCov = 
  jacobianLocalToGlobal.transpose() * rotatedMatrix * jacobianLocalToGlobal;

  printMatrix("Rotated to Acts local cov : ",actsLocalCov);
  return actsLocalCov;

}


Acts::BoundSymMatrix ActsTransformations::rotateActsCovToSvtxTrack(
			        const Acts::BoundTrackParameters params,
				Acts::GeometryContext geoCtxt) const
{

  auto covarianceMatrix = *params.covariance();

  printMatrix("Initial Acts covariance: ", covarianceMatrix);

  const double px = params.momentum()(0);
  const double py = params.momentum()(1);
  const double pz = params.momentum()(2);
  const double p = params.momentum().norm();
  
  const double uPx = px / p;
  const double uPy = py / p;
  const double uPz = pz / p;
  
  const double cosTheta = uPz;
  const double sinTheta = sqrt(uPx * uPx + uPy * uPy);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = uPx * invSinTheta;
  const double sinPhi = uPy * invSinTheta;

  Acts::BoundToFreeMatrix jacobianLocalToGlobal = Acts::BoundToFreeMatrix::Zero();
  jacobianLocalToGlobal(0, Acts::eBoundLoc0) = -sinPhi;
  jacobianLocalToGlobal(0, Acts::eBoundLoc1) = -cosPhi * cosTheta;
  jacobianLocalToGlobal(1, Acts::eBoundLoc0) = cosPhi;
  jacobianLocalToGlobal(1, Acts::eBoundLoc1) = -sinPhi * cosTheta;
  jacobianLocalToGlobal(2, Acts::eBoundLoc1) = sinTheta;
  jacobianLocalToGlobal(3, Acts::eBoundTime) = 1;
  jacobianLocalToGlobal(4, Acts::eBoundPhi) = -sinTheta * sinPhi;
  jacobianLocalToGlobal(4, Acts::eBoundTheta) = cosTheta * cosPhi;
  jacobianLocalToGlobal(5, Acts::eBoundPhi) = sinTheta * cosPhi;
  jacobianLocalToGlobal(5, Acts::eBoundTheta) = cosTheta * sinPhi;
  jacobianLocalToGlobal(6, Acts::eBoundTheta) = -sinTheta;
  jacobianLocalToGlobal(7, Acts::eBoundQOverP) = 1;

  /// Covariance is now an 8x8 matrix in basis (x,y,z,time,Tx,Ty,Tz,q/p)
  auto rotatedMatrix 
    = jacobianLocalToGlobal * covarianceMatrix * jacobianLocalToGlobal.transpose();

  /// Now rotate to x,y,z, px,py,pz
  /// ActsMatrixD is an eigen matrix
  Acts::ActsMatrixD<6,8> sphenixRot;
  sphenixRot.setZero();
  
  /// Make the xyz transform unity
  sphenixRot(0,0) = 1;
  sphenixRot(1,1) = 1;
  sphenixRot(2,2) = 1;
  sphenixRot(3,4) = p;
  sphenixRot(4,5) = p;
  sphenixRot(5,6) = p;
  sphenixRot(3,7) = uPx * p * p;
  sphenixRot(4,7) = uPy * p * p;
  sphenixRot(5,7) = uPz * p * p;
  
  Acts::BoundSymMatrix globalCov = Acts::BoundSymMatrix::Zero();
  globalCov = sphenixRot * rotatedMatrix * sphenixRot.transpose();

  printMatrix("Global sPHENIX cov : ", globalCov);

  /// Convert to sPHENIX units
  for(int i = 0; i < 6; ++i)
    {
      for(int j = 0; j < 6; ++j)
	{
	  if(i < 3 && j < 3)
	    globalCov(i,j) /= Acts::UnitConstants::cm2;
	  else if (i < 3)
	    globalCov(i,j) /= Acts::UnitConstants::cm;
	  else if (j < 3)
	    globalCov(i,j) /= Acts::UnitConstants::cm;
	  
	}
    }

  printMatrix("Global sphenix cov after unit conv: " , globalCov);

  return globalCov;
}


void ActsTransformations::printMatrix(const std::string &message,
					Acts::BoundSymMatrix matrix) const
{
 
  if(m_verbosity > 10)
    {
      std::cout << std::endl << message.c_str() << std::endl;
      for(int i = 0 ; i < matrix.rows(); ++i)
	{
	  std::cout << std::endl;
	  for(int j = 0; j < matrix.cols(); ++j)
	    {
	      std::cout << matrix(i,j) << ", ";
	    }
	}
    
      std::cout << std::endl;
    }
}



void ActsTransformations::calculateDCA(const Acts::BoundTrackParameters param,
				       Acts::Vector3D vertex,
				       Acts::BoundSymMatrix cov,
				       Acts::GeometryContext geoCtxt,
				       float &dca3Dxy,
				       float &dca3Dz,
				       float &dca3DxyCov,
				       float &dca3DzCov) const
{
  Acts::Vector3D pos = param.position(geoCtxt);
  Acts::Vector3D mom = param.momentum();

  /// Correct for initial vertex estimation
  pos -= vertex;

  Acts::ActsSymMatrixD<3> posCov;
  for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
	{
	  posCov(i,j) = cov(i,j);
	} 
    }

  Acts::Vector3D r = mom.cross(Acts::Vector3D(0.,0.,1.));
  float phi = atan2(r(1), r(0));

  Acts::RotationMatrix3D rot;
  Acts::RotationMatrix3D rot_T;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  
  rot_T = rot.transpose();

  Acts::Vector3D pos_R = rot * pos;
  Acts::ActsSymMatrixD<3> rotCov = rot * posCov * rot_T;

  dca3Dxy = pos_R(0);
  dca3Dz = pos_R(2);
  dca3DxyCov = rotCov(0,0);
  dca3DzCov = rotCov(2,2);
  
}



void ActsTransformations::fillSvtxTrackStates(const Trajectory& traj,
					      const size_t &trackTip,
					      SvtxTrack *svtxTrack,
					      Acts::GeometryContext geoContext) const
{

 const auto &[trackTips, mj] = traj.trajectory();
  
  mj.visitBackwards(trackTip, [&](const auto &state) 
  {
    
      /// Only fill the track states with non-outlier measurement
      const auto typeFlags = state.typeFlags();
      if( !typeFlags.test(Acts::TrackStateFlag::MeasurementFlag) )
      { return true; }
      
      // only fill for state vectors with proper smoothed parameters
      if( !state.hasSmoothed()) return true;

      // create svtx state vector with relevant pathlength
      const float pathlength = state.pathLength() / Acts::UnitConstants::cm;  
      SvtxTrackState_v1 out( pathlength );
    
      // get smoothed fitted parameters
      const Acts::BoundTrackParameters params(state.referenceSurface().getSharedPtr(),
        state.smoothed(),
        state.smoothedCovariance());
      
      // position
      const auto global = params.position(geoContext);
      out.set_x(global.x() / Acts::UnitConstants::cm);
      out.set_y(global.y() / Acts::UnitConstants::cm);
      out.set_z(global.z() / Acts::UnitConstants::cm);
      
      // momentum
      const auto momentum = params.momentum();
      out.set_px(momentum.x());
      out.set_py(momentum.y());
      out.set_pz(momentum.z());
      
      /// covariance    
      const auto globalCov = rotateActsCovToSvtxTrack(params, geoContext);
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
      { out.set_error(i, j, globalCov(i,j)); }

      // add cluster key
      const auto cluskey = state.uncalibrated().cluskey();
      svtxTrack->insert_cluster_key(cluskey);

      // print
      if(m_verbosity > 20)
      {
        std::cout << " inserting state with x,y,z ="
          << " " << global.x() /  Acts::UnitConstants::cm 
          << " " << global.y() /  Acts::UnitConstants::cm 
          << " " << global.z() /  Acts::UnitConstants::cm 
          << " pathlength " << pathlength
          << " momentum px,py,pz = " <<  momentum.x() << "  " <<  momentum.y() << "  " << momentum.y()  
          << " cluskey " << cluskey << std::endl
          << "covariance " << globalCov << std::endl; 
      }
	  
      svtxTrack->insert_state(&out);      
  
      return true;      
    }
    );

  return;
}
