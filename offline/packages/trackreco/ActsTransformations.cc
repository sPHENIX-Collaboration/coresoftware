#include "ActsTransformations.h"
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <chrono>
using namespace std::chrono;

Acts::BoundSymMatrix ActsTransformations::rotateSvtxTrackCovToActs(
			        const SvtxTrack *track,
				Acts::GeometryContext geoCtxt)
{
  Acts::BoundSymMatrix rotation = Acts::BoundSymMatrix::Zero();
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
  const int charge = track->get_charge();
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
    
  const double x = track->get_x();
  const double y = track->get_y();
  const double z = track->get_z();
  const double r = sqrt(x*x + y*y + z*z);
  
  const double posCosTheta = z / r;
  const double posSinTheta = sqrt(x*x + y*y) / r;
  const double posInvSinTheta = 1. / posSinTheta;
  const double posCosPhi = x * posInvSinTheta / r;
  const double posSinPhi = y * posInvSinTheta / r;

  /// Position rotation to Acts loc0 and loc1, which are the local points
  /// on a surface centered at the (x,y,z) global position with normal
  /// vector in the direction of the unit momentum vector
  rotation(0,0) = - posSinPhi;
  rotation(0,1) =   posCosPhi;
  rotation(1,0) =   posCosPhi * posCosTheta;
  rotation(1,1) = - posSinPhi * posCosTheta;
  rotation(1,2) = - posSinTheta;

  // Directional and momentum parameters for curvilinear
  rotation(2, 3) = - p * sinPhi * sinTheta;
  rotation(2, 4) =   p * cosPhi * sinTheta;
  rotation(3, 3) =   p * cosPhi * cosTheta;
  rotation(3, 4) =   p * sinPhi * cosTheta;
  rotation(3, 5) = - p * sinTheta;

  ///q/p rotaton
  rotation(4,3) = -charge * px / pow(p,0.5);
  rotation(4,4) = -charge * py / pow(p,0.5);
  rotation(4,5) = -charge * pz / pow(p,0.5);

  /// time is left as 0
  printMatrix("rotation is : ",rotation);
  Acts::BoundSymMatrix matrix = rotation * svtxCovariance * rotation.transpose();
  printMatrix("Rotated is : ",matrix);
  return matrix;

}


Acts::BoundSymMatrix ActsTransformations::rotateActsCovToSvtxTrack(
			        const Acts::BoundTrackParameters params,
				Acts::GeometryContext geoCtxt)
{

  auto covarianceMatrix = *params.covariance();
  
  printMatrix("Initial Acts covariance: ", covarianceMatrix);

  const double px = params.momentum()(0);
  const double py = params.momentum()(1);
  const double pz = params.momentum()(2);
  const double p = sqrt(px * px + py * py + pz * pz);
  
  const double x = params.position(geoCtxt)(0);
  const double y = params.position(geoCtxt)(1);
  const double z = params.position(geoCtxt)(2);
  const double r = sqrt(x*x + y*y + z*z);

  const double posCosTheta = z / r;
  const double posSinTheta = sqrt(x*x + y*y) / r;
  const double posCosPhi = x / ( r * posSinTheta);
  const double posSinPhi = y / ( r * posSinTheta);
  
  const int charge = params.charge();
  
  const double uPx = px / p;
  const double uPy = py / p;
  const double uPz = pz / p;
  
  const double cosTheta = uPz;
  const double sinTheta = sqrt(uPx * uPx + uPy * uPy);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = uPx * invSinTheta;
  const double sinPhi = uPy * invSinTheta;

  Acts::BoundSymMatrix rotation = Acts::BoundSymMatrix::Zero();

  /// This is the original matrix we rotated by. So instead of rotating
  /// as normal RCR^T, we will just rotate back by rotating by the same 
  /// matrix as R^TCR
  /// Position rotation to Acts loc0 and loc1, which are the local points
  /// on a surface centered at the (x,y,z) global position with normal
  /// vector in the direction of the unit momentum vector
  rotation(0,0) = - posSinPhi;
  rotation(0,1) =   posCosPhi;
  rotation(1,0) = - posCosPhi * posCosTheta;
  rotation(1,1) = - posSinPhi * posCosTheta;
  rotation(1,2) =   posSinTheta;

  // Directional and momentum parameters for curvilinear
  rotation(2, 3) = -p * sinPhi * sinTheta;
  rotation(2, 4) =  p * cosPhi * sinTheta;
  rotation(3, 3) =  p * cosPhi * cosTheta;
  rotation(3, 4) =  p * sinPhi * cosTheta;
  rotation(3, 5) = -p * sinTheta;
  
  ///q/p rotaton
  ///d(q/p)/dp_i = q * -p_i * p^{-3/2}
  rotation(4,3) = -charge * px / pow(p,1.5);
  rotation(4,4) = -charge * py / pow(p,1.5);
  rotation(4,5) = -charge * pz / pow(p,1.5);

  printMatrix("Rotating back to global with : ", rotation.transpose());

  Acts::BoundSymMatrix globalCov = Acts::BoundSymMatrix::Zero();
  globalCov = rotation.transpose() * covarianceMatrix * rotation;


  /// convert back to sPHENIX coordinates of cm
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

  printMatrix("Global sPHENIX cov : ", globalCov);


  return globalCov;
}


void ActsTransformations::printMatrix(const std::string &message,
					Acts::BoundSymMatrix matrix)
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
				   Acts::GeometryContext geoCtxt,
				   float &dca3Dxy,
				   float &dca3Dz,
				   float &dca3DxyCov,
				   float &dca3DzCov)
{
  Acts::Vector3D pos = param.position(geoCtxt);
  Acts::Vector3D mom = param.momentum();

  /// Correct for initial vertex estimation
  pos -= vertex;

  Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
  if(param.covariance())
    cov = param.covariance().value();

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



void ActsTransformations::fillSvtxTrackStates(const Trajectory traj,
					      const size_t &trackTip,
					      SvtxTrack *svtxTrack,
					      Acts::GeometryContext geoContext,
					      CluskeyBimap *hitIDCluskeyMap)
{

 const auto &[trackTips, mj] = traj.trajectory();
  
  mj.visitBackwards(trackTip, [&](const auto &state) {
      /// Only fill the track states with non-outlier measurement
      auto typeFlags = state.typeFlags();
      if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag))
	{
	  return true;
	}
      
      auto meas = std::get<Measurement>(*state.uncalibrated());

      /// Get local position
      Acts::Vector2D local(meas.parameters()[Acts::eBoundLoc0],
			   meas.parameters()[Acts::eBoundLoc1]);

      /// This is an arbitrary vector. Doesn't matter in coordinate transformation
      /// in Acts code
      Acts::Vector3D mom(1., 1., 1.);
      Acts::Vector3D global = meas.referenceObject().localToGlobal(
					    geoContext,
					    local, mom);
      
      float pathlength = state.pathLength() / Acts::UnitConstants::cm;  
      SvtxTrackState_v1 out( pathlength );
      out.set_x(global.x() / Acts::UnitConstants::cm);
      out.set_y(global.y() / Acts::UnitConstants::cm);
      out.set_z(global.z() / Acts::UnitConstants::cm);
    
      if (state.hasSmoothed())
	{
	  Acts::BoundTrackParameters parameter(state.referenceSurface().getSharedPtr(),
					  state.smoothed(),
					  state.smoothedCovariance());
	  out.set_px(parameter.momentum().x());
	  out.set_py(parameter.momentum().y());
	  out.set_pz(parameter.momentum().z());

	  /// Get measurement covariance    

	  Acts::BoundSymMatrix globalCov = rotateActsCovToSvtxTrack(parameter,
								    geoContext);
	  for (int i = 0; i < 6; i++)
	    {
	      for (int j = 0; j < 6; j++)
		{ 
		  out.set_error(i, j, globalCov(i,j)); 
		}
	    }

	  const unsigned int hitId = state.uncalibrated().hitID();
	  TrkrDefs::cluskey cluskey = 
	    hitIDCluskeyMap->right.find(hitId)->first;
	  svtxTrack->insert_cluster_key(cluskey);

	  if(m_verbosity > 20)
	    {
	      std::cout << " inserting state with x,y,z = " << global.x() /  Acts::UnitConstants::cm 
			<< "  " << global.y() /  Acts::UnitConstants::cm << "  " 
			<< global.z() /  Acts::UnitConstants::cm 
			<< " pathlength " << pathlength
			<< " momentum px,py,pz = " <<  parameter.momentum().x() << "  " <<  parameter.momentum().y() << "  " << parameter.momentum().y()  
			<< " cluskey " << cluskey << std::endl
			<< "covariance " << globalCov << std::endl; 
	    }
	  
	  svtxTrack->insert_state(&out);      
	}
  
      return true;      
    }
    );

  return;
}
