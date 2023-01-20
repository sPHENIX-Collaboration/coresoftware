#include "ActsTransformations.h"
#include "SvtxTrack.h"
#include "SvtxTrackState.h"
#include "SvtxTrackState_v1.h"

#include <trackbase/TrkrCluster.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>

namespace
{
  /// square
  template<class T> inline constexpr T square(const T& x) {return x*x;}

  /// get radius from coordinates
  template<class T> T radius(const T& x, const T& y)
  { return std::sqrt(square(x) + square(y));}

  /// templatize print matrix method
  template<class T>
    void print_matrix( const std::string &message, const T& matrix)
  {
    std::cout << std::endl << message << std::endl;
    for(int i = 0 ; i < matrix.rows(); ++i)
    {
      for(int j = 0; j < matrix.cols(); ++j)
      { std::cout << matrix(i,j) << ", "; }
      
      std::cout << std::endl;
    }
  }
    
}

//_______________________________________________________________________________
Acts::BoundSymMatrix ActsTransformations::rotateSvtxTrackCovToActs( const SvtxTrack *track ) const
{ return rotateSvtxTrackCovToActs( track->find_state(0.0)->second ); }

//_______________________________________________________________________________
Acts::BoundSymMatrix ActsTransformations::rotateSvtxTrackCovToActs(
			        const SvtxTrackState *state) const
{
  Acts::BoundSymMatrix svtxCovariance = Acts::BoundSymMatrix::Zero();

  for(int i = 0; i < 6; ++i)
    {
      for(int j = 0; j < 6; ++j)
	{
	  svtxCovariance(i,j) = state->get_error(i,j);
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
 
  const double px = state->get_px();
  const double py = state->get_py();
  const double pz = state->get_pz();
  const double p2 = square(px) + square(py) + square(pz);
  const double p = std::sqrt(p2);
  const double invp = 1./p;
  
  const double uPx = px/p;
  const double uPy = py/p;
  const double uPz = pz/p;

  //Acts version
  const double cosTheta = uPz;
  const double sinTheta = std::sqrt(square(uPx) + square(uPy));
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = uPx * invSinTheta; // equivalent to x/r
  const double sinPhi = uPy * invSinTheta; // equivalent to y/r
    
  /// First we rotate to (x,y,z,time,Tx,Ty,Tz,q/p) to take advantage of the
  /// already created Acts rotation matrix from this basis into the Acts local basis
  /// We basically go backwards from rotateActsCovToSvtxTrack to get the Acts cov 
  /// from the SvtxTrack cov

  // rotate to (x,y,z,t,upx,upy,upz,q/p)
  Acts::ActsMatrix<8,6> sphenixRot = Acts::ActsMatrix<8,6>::Zero(); 
  sphenixRot(0,0) = 1;
  sphenixRot(1,1) = 1;
  sphenixRot(2,2) = 1;
  sphenixRot(4,3) = invp;
  sphenixRot(5,4) = invp;
  sphenixRot(6,5) = invp;
  sphenixRot(7,3) = uPx/p2;
  sphenixRot(7,4) = uPy/p2;
  sphenixRot(7,5) = uPz/p2;

  const auto rotatedMatrix = sphenixRot*svtxCovariance*sphenixRot.transpose();
  
  // global to local transformation
  /* from https://github.com/acts-project/acts/blob/394cfdb308956de93f90ab3162040bb6d835027d/Core/src/Propagator/detail/CovarianceEngine.cpp#L31 */
  Acts::FreeToBoundMatrix jacobianGlobalToLocal = Acts::FreeToBoundMatrix::Zero();
  jacobianGlobalToLocal(Acts::eBoundLoc0, 0) = -sinPhi;
  jacobianGlobalToLocal(Acts::eBoundLoc0, 1) = cosPhi;
  jacobianGlobalToLocal(Acts::eBoundLoc1, 0) = -cosPhi * cosTheta;
  jacobianGlobalToLocal(Acts::eBoundLoc1, 1) = -sinPhi * cosTheta;
  jacobianGlobalToLocal(Acts::eBoundLoc1, 2) = sinTheta;
  jacobianGlobalToLocal(Acts::eBoundTime, 3) = 1.;
  jacobianGlobalToLocal(Acts::eBoundPhi, 4) = -sinPhi * invSinTheta;
  jacobianGlobalToLocal(Acts::eBoundPhi, 5) = cosPhi * invSinTheta;
  jacobianGlobalToLocal(Acts::eBoundTheta, 4) = cosPhi * cosTheta;
  jacobianGlobalToLocal(Acts::eBoundTheta, 5) = sinPhi * cosTheta;
  jacobianGlobalToLocal(Acts::eBoundTheta, 6) = -sinTheta;
  jacobianGlobalToLocal(Acts::eBoundQOverP, 7) = 1.;
  
  Acts::BoundSymMatrix actsLocalCov = jacobianGlobalToLocal * rotatedMatrix * jacobianGlobalToLocal.transpose();

  /*
   * need to assign the covariance matrix diagonal element corresponding to the time coordinate manually, 
   * since it is lost in the conversion to Svtx coordinates 
   */
  actsLocalCov( Acts::eBoundTime, Acts::eBoundTime ) = 1.;
  
  printMatrix("Rotated to Acts local cov : ",actsLocalCov);
  return actsLocalCov;

}

//_______________________________________________________________________________
Acts::BoundSymMatrix ActsTransformations::rotateActsCovToSvtxTrack( const ActsTrackFittingAlgorithm::TrackParameters& params ) const
{ 
  const auto covarianceMatrix = *params.covariance();
  printMatrix("Initial Acts covariance: ", covarianceMatrix);

  const double px = params.momentum().x();
  const double py = params.momentum().y();
  const double pz = params.momentum().z();
  const double p = params.momentum().norm();
  const double p2 = square(p);
  
  const double uPx = px/p;
  const double uPy = py/p;
  const double uPz = pz/p;
  
  const double cosTheta = uPz;
  const double sinTheta = std::sqrt(square(uPx) + square(uPy));
  const double invSinTheta = 1./sinTheta;
  const double cosPhi = uPx*invSinTheta;
  const double sinPhi = uPy*invSinTheta;

  // local to global transformation
  /* from https://github.com/acts-project/acts/blob/394cfdb308956de93f90ab3162040bb6d835027d/Core/src/Propagator/detail/CovarianceEngine.cpp#L222 */
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
  
  // Covariance is now an 8x8 matrix in basis (x,y,z,time,Tx,Ty,Tz,q/p)
  const auto rotatedMatrix = jacobianLocalToGlobal * covarianceMatrix * jacobianLocalToGlobal.transpose();

  // Now rotate to x,y,z, px,py,pz
  Acts::ActsMatrix<6,8> sphenixRot = Acts::ActsMatrix<6,8>::Zero(); 
  sphenixRot(0,0) = 1;
  sphenixRot(1,1) = 1;
  sphenixRot(2,2) = 1;
  sphenixRot(3,4) = p;
  sphenixRot(4,5) = p;
  sphenixRot(5,6) = p;
  sphenixRot(3,7) = uPx*p2;
  sphenixRot(4,7) = uPy*p2;
  sphenixRot(5,7) = uPz*p2;
  
  Acts::BoundSymMatrix globalCov = sphenixRot * rotatedMatrix * sphenixRot.transpose();
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

void ActsTransformations::printMatrix(const std::string &message, const Acts::BoundSymMatrix& matrix) const
{ if(m_verbosity > 10) print_matrix( message, matrix ); }

void ActsTransformations::calculateDCA(const Acts::BoundTrackParameters param,
				       Acts::Vector3 vertex,
				       Acts::BoundSymMatrix cov,
				       Acts::GeometryContext& geoCtxt,
				       float &dca3Dxy,
				       float &dca3Dz,
				       float &dca3DxyCov,
				       float &dca3DzCov) const
{
  Acts::Vector3 pos = param.position(geoCtxt);
  Acts::Vector3 mom = param.momentum();

  /// Correct for initial vertex estimation
  pos -= vertex;

  Acts::ActsSymMatrix<3> posCov;
  for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
	{
	  posCov(i,j) = cov(i,j);
	} 
    }

  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));
  float phi = atan2(r(1), r(0));

  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
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

  Acts::Vector3 pos_R = rot * pos;
  Acts::ActsSymMatrix<3> rotCov = rot * posCov * rot_T;

  dca3Dxy = pos_R(0);
  dca3Dz = pos_R(2);
  dca3DxyCov = rotCov(0,0);
  dca3DzCov = rotCov(2,2);
  
}


void ActsTransformations::fillSvtxTrackStates(const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>& traj,
					      const size_t& trackTip,
					      SvtxTrack *svtxTrack,
					      Acts::GeometryContext& geoContext) const
{ 
  traj.visitBackwards(trackTip, [&](const auto &state) 
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
      const auto globalCov = rotateActsCovToSvtxTrack(params);
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
      { out.set_error(i, j, globalCov(i,j)); }

      // print
      if(m_verbosity > 20)
      {
        std::cout << " inserting state with x,y,z ="
          << " " << global.x() /  Acts::UnitConstants::cm 
          << " " << global.y() /  Acts::UnitConstants::cm 
          << " " << global.z() /  Acts::UnitConstants::cm 
          << " pathlength " << pathlength
          << " momentum px,py,pz = " <<  momentum.x() << "  " <<  momentum.y() << "  " << momentum.y()  
		  << std::endl
          << "covariance " << globalCov << std::endl; 
      }

      svtxTrack->insert_state(&out);      
  
      return true;      
    }
    );

  return;
}
