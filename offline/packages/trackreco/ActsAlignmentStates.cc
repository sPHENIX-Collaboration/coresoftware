#include "ActsAlignmentStates.h"
#include "ActsPropagator.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSourceLink.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxAlignmentStateMap.h>
#include <trackbase_historic/SvtxAlignmentState_v1.h>
#include <trackbase_historic/SvtxTrack.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/Surfaces/Surface.hpp>

namespace
{
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

void ActsAlignmentStates::fillAlignmentStateMap(const Trajectory& traj,
                                                SvtxTrack* track,
						const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  const auto mj = traj.multiTrajectory();
  const auto& tips = traj.tips();
  const auto& trackTip = tips.front();
  const auto crossing = track->get_silicon_seed()->get_crossing();

  /// trajectory (helix) paramters
  const auto& params = traj.trackParameters(trackTip);
  const auto boundparams = params.parameters();

  ActsPropagator propagator(m_tGeometry);

  SvtxAlignmentStateMap::StateVec statevec;
  if(m_verbosity > 2) std::cout << "Beginning alignment state creation for track " << track->get_id() << std::endl;
  std::vector<Acts::BoundIndices> indices {Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi, Acts::eBoundTheta, Acts::eBoundQOverP, Acts::eBoundTime};
  mj.visitBackwards(trackTip, [&](const auto& state)
                    {
    /// Collect only track states which were used in smoothing of KF and are measurements
    if (not state.hasSmoothed() or
        not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag))
    {
      return true;
    }
      
    const auto& surface = state.referenceSurface();
    const auto& sl = static_cast<const ActsSourceLink&>(state.uncalibrated());
    auto ckey = sl.cluskey();
    Acts::Vector2 localMeas = Acts::Vector2::Zero();
    /// get the local measurement that acts used
    std::visit([&](const auto& meas) {
	localMeas(0) = meas.parameters()[0];
	localMeas(1) = meas.parameters()[1];
      }, measurements[sl.index()]);
    
    if (m_verbosity > 2)
      {
	std::cout << "sl index and ckey " << sl.index() << ", "
		  << sl.cluskey() << " with local position " 
		  << localMeas.transpose() << std::endl;
      }

    auto clus = m_clusterMap->findCluster(ckey);
    const auto trkrId = TrkrDefs::getTrkrId(ckey);
 
    const Acts::Vector2 localState = state.effectiveProjector() * state.smoothed();
    /// Local residual between measurement and smoothed Acts state
    const Acts::Vector2 localResidual = localMeas - localState;  

    Acts::Vector3 clusGlobal = m_tGeometry->getGlobalPosition(ckey, clus);
    if (trkrId == TrkrDefs::tpcId)
    {
      makeTpcGlobalCorrections(ckey, crossing, clusGlobal);
    }
   
    /// convert to acts units
    clusGlobal *= Acts::UnitConstants::cm;

    const Acts::FreeVector globalStateParams = Acts::detail::transformBoundToFreeParameters(surface, m_tGeometry->geometry().getGeoContext(), state.smoothed());
    Acts::Vector3 stateGlobal = globalStateParams.segment<3>(Acts::eFreePos0);

    Acts::Vector3 clus_sigma(0, 0, 0);

    if (m_clusterVersion != 4)
      {
        clus_sigma(2) = clus->getZError() * Acts::UnitConstants::cm;
        clus_sigma(0) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
        clus_sigma(1) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
      }
    else 
      {
        double clusRadius = sqrt(clusGlobal(0) * clusGlobal(0) + clusGlobal(1) * clusGlobal(1)) / Acts::UnitConstants::cm;
        auto para_errors = m_clusErrPara.get_simple_cluster_error(clus, clusRadius, ckey);
        float exy2 = para_errors.first * Acts::UnitConstants::cm2;
        float ez2 = para_errors.second * Acts::UnitConstants::cm2;
        clus_sigma(2) = sqrt(ez2);
        clus_sigma(0) = sqrt(exy2 / 2.0);
        clus_sigma(1) = sqrt(exy2 / 2.0);
      }
  
 
    if (m_verbosity > 2)
    {
      std::cout << "clus global is " << clusGlobal.transpose() << std::endl
                << "state global is " << stateGlobal.transpose() << std::endl;
      std::cout << "   clus errors are " << clus_sigma.transpose() << std::endl;
      std::cout << "   clus loc is " << localMeas.transpose() << std::endl;
      std::cout << "   state loc is " << localState << std::endl;
      std::cout << "   local residual is " << localResidual.transpose() << std::endl;
    }

    // Get the derivative of alignment (global) parameters w.r.t. measurement or residual
    /// The local bound parameters still have access to global phi/theta
    double phi = state.smoothed()[Acts::eBoundPhi];
    double theta = state.smoothed()[Acts::eBoundTheta];
    const Acts::Vector3 tangent = Acts::makeDirectionUnitFromPhiTheta(phi,theta);
    if(m_verbosity > 2)
      {
	std::cout << "tangent vector to track state is " << tangent.transpose() << std::endl;
      }

    std::pair<Acts::Vector3, Acts::Vector3> projxy = 
      get_projectionXY(surface, tangent);

    Acts::Vector3 sensorCenter = surface.center(m_tGeometry->geometry().getGeoContext());
    Acts::Vector3 OM = stateGlobal - sensorCenter;
    if(m_verbosity > 2)
      {
	std::cout << "   global deriv calcs" << std::endl
		  << "stateGlobal: " << stateGlobal.transpose()
		  << ", sensor center " << sensorCenter.transpose() << std::endl
		  << ", OM " << OM.transpose() << std::endl << "   projxy "
		  << projxy.first.transpose() << ", " 
		  << projxy.second.transpose() << std::endl;
      }

    //! this is the derivative of the state wrt to Acts track parameters
    //! e.g. (d_0, z_0, phi, theta, q/p, t)
  
    auto localDeriv = -state.effectiveProjector() * state.jacobian();
    if(m_verbosity > 2)
      {
	std::cout << "local deriv " << std::endl << localDeriv << std::endl;
	std::cout << "local deriv rows cols " << localDeriv.rows() << ", " << localDeriv.cols() << std::endl;
      }

    SvtxAlignmentState::LocalMatrix locDerivCalc = SvtxAlignmentState::LocalMatrix::Zero();
    for(int nloc = 0; nloc < SvtxAlignmentState::NLOC; nloc++)
      {
	Acts::BoundVector newparam;
	for(auto index : indices)
	  {
	    newparam[index] = boundparams[index];
	    if(index == nloc) newparam[index ] = boundparams[index]+0.1;
	  }
	const Acts::BoundTrackParameters newparams(
             params.referenceSurface().getSharedPtr(), 
	     newparam, params.charge(), params.covariance());
	propagator.verbosity(m_verbosity);
	/// propagate the helix to the surface with the modified track parameters
	auto result = propagator.propagateTrack(newparams, surface.getSharedPtr());
	if(result.ok())
	  {
	    const Acts::Vector3 newStateGlob = result.value().second.position(m_tGeometry->geometry().getGeoContext());
	    if(m_verbosity > 2)
	      {
		std::cout << "  kf state " << stateGlobal.transpose()
			  << std::endl << "  new state " << newStateGlob.transpose() << std::endl;
	      }
	    Acts::Vector3 stateDiff = newStateGlob - stateGlobal;
	    stateDiff /= 0.1;
	    locDerivCalc(0, nloc) = stateDiff.dot(projxy.first);
	    locDerivCalc(1, nloc) = stateDiff.dot(projxy.second);
	    
	  }
      }
    if(m_verbosity > 2)
      {
	std::cout << "hand calculated local derivatives " 
		  << std::endl << locDerivCalc << std::endl;

      }

    auto globDeriv = makeGlobalDerivatives(OM, projxy);
    auto svtxstate = std::make_unique<SvtxAlignmentState_v1>();
    
    svtxstate->set_residual(localResidual);
    svtxstate->set_local_derivative_matrix(localDeriv);
    svtxstate->set_global_derivative_matrix(globDeriv);
   svtxstate->set_cluster_key(ckey);
  
    statevec.push_back(svtxstate.release());
   return true; });

  if (m_verbosity > 2)
  {
    std::cout << "Inserting track " << track->get_id() << " with nstates "
              << statevec.size() << std::endl;
  }

  m_alignmentStateMap->insertWithKey(track->get_id(), statevec);

  return;
}

SvtxAlignmentState::GlobalMatrix 
ActsAlignmentStates::makeGlobalDerivatives(const Acts::Vector3& OM, 
					   const std::pair<Acts::Vector3,Acts::Vector3>& projxy)
{
  SvtxAlignmentState::GlobalMatrix globalder = SvtxAlignmentState::GlobalMatrix::Zero();
  Acts::SymMatrix3 unit = Acts::SymMatrix3::Identity();

  //! x residual rotations
  globalder(0,0) = ((unit.col(0)).cross(OM)).dot(projxy.first);
  globalder(0,1) = ((unit.col(1)).cross(OM)).dot(projxy.first);
  globalder(0,2) = ((unit.col(2)).cross(OM)).dot(projxy.first);
  //! x residual translations
  globalder(0,3) = unit.col(0).dot(projxy.first);
  globalder(0,4) = unit.col(1).dot(projxy.first);
  globalder(0,5) = unit.col(2).dot(projxy.first);
  
  //! y residual rotations
  globalder(1,0) = ((unit.col(0)).cross(OM)).dot(projxy.second);
  globalder(1,1) = ((unit.col(1)).cross(OM)).dot(projxy.second);
  globalder(1,2) = ((unit.col(2)).cross(OM)).dot(projxy.second);
  //! y residual translations
  globalder(1,3) = unit.col(0).dot(projxy.second);
  globalder(1,4) = unit.col(1).dot(projxy.second);
  globalder(1,5) = unit.col(2).dot(projxy.second);
  
  return globalder;

}

std::pair<Acts::Vector3, Acts::Vector3> ActsAlignmentStates::get_projectionXY(const Acts::Surface& surface, const Acts::Vector3& tangent)
{
  Acts::Vector3 projx = Acts::Vector3::Zero();
  Acts::Vector3 projy = Acts::Vector3::Zero();
  
  // get the plane of the surface
  Acts::Vector3 sensorCenter = surface.center(m_tGeometry->geometry().getGeoContext()); 
  // sensorNormal is the Z vector
  Acts::Vector3 Z = -surface.normal(m_tGeometry->geometry().getGeoContext());

  // get surface X and Y unit vectors in global frame
  // transform Xlocal = 1.0 to global, subtract the surface center, normalize to 1
  Acts::Vector3 xloc(1.0,0.0,0.0);
  Acts::Vector3 xglob =  surface.transform(m_tGeometry->geometry().getGeoContext()) * xloc;

  Acts::Vector3 yloc(0.0,1.0,0.0);
  Acts::Vector3 yglob =  surface.transform(m_tGeometry->geometry().getGeoContext()) * yloc;

  Acts::Vector3 X = (xglob-sensorCenter) / (xglob-sensorCenter).norm();
  Acts::Vector3 Y = (yglob-sensorCenter) / (yglob-sensorCenter).norm();

  projx = X - (tangent.dot(X) / tangent.dot(Z)) * Z;
  projy = Y - (tangent.dot(Y) / tangent.dot(Z)) * Z;

  return std::make_pair(projx, projy);
}

void ActsAlignmentStates::getGlobalDerivatives(std::vector<Acts::Vector3>& anglederivs, SvtxAlignmentState::GlobalMatrix& analytic)
{
  analytic(0,0) = 1;
  analytic(1,1) = 1;
  analytic(2,2) = 1;
  
  for(int res = 0; res < SvtxAlignmentState::NRES; ++res)
    {
      for(int gl = 3; gl < SvtxAlignmentState::NGL; ++gl)
	{
	  // this is not backwards - the angle derivs is transposed
	  analytic(res,gl) = anglederivs.at(gl-3)(res) * Acts::UnitConstants::cm;
	}
    }
  
}

void ActsAlignmentStates::makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global)
{
  // make all corrections to global position of TPC cluster
  unsigned int side = TpcDefs::getSide(cluster_key);
  float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
  global[2] = z;

  // apply distortion corrections
  if (m_dcc_static)
  {
    global = m_distortionCorrection.get_corrected_position(global, m_dcc_static);
  }
  if (m_dcc_average)
  {
    global = m_distortionCorrection.get_corrected_position(global, m_dcc_average);
  }
  if (m_dcc_fluctuation)
  {
    global = m_distortionCorrection.get_corrected_position(global, m_dcc_fluctuation);
  }
}

std::vector<Acts::Vector3> ActsAlignmentStates::getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing)
{
  // The value of global is from the geocontext transformation
  // we add to that transformation a small rotation around the relevant axis in the surface frame

  std::vector<Acts::Vector3> derivs_vector;

  // get the transformation from the geocontext
  Acts::Transform3 transform = surface->transform(m_tGeometry->geometry().getGeoContext());

  // Make an additional transform that applies a small rotation angle in the surface frame
  for (unsigned int iangle = 0; iangle < 3; ++iangle)
  {
    // creates transform that adds a perturbation rotation around one axis, uses it to estimate derivative wrt perturbation rotation

    unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
    unsigned int layer = TrkrDefs::getLayer(cluster_key);

    Acts::Vector3 derivs(0, 0, 0);
    Eigen::Vector3d theseAngles(0, 0, 0);
    theseAngles[iangle] = sensorAngles[iangle];  // set the one we want to be non-zero

    Acts::Vector3 keeper(0, 0, 0);
    for (int ip = 0; ip < 2; ++ip)
    {
      if (ip == 1)
      {
        theseAngles[iangle] *= -1;
      }  // test both sides of zero

      if (m_verbosity > 1)
      {
        std::cout << "     trkrId " << trkrId << " layer " << layer << " cluster_key " << cluster_key
                  << " sensorAngles " << theseAngles[0] << "  " << theseAngles[1] << "  " << theseAngles[2] << std::endl;
      }

      Acts::Transform3 perturbationTransformation = makePerturbationTransformation(theseAngles);
      Acts::Transform3 overallTransformation = transform * perturbationTransformation;

      // transform the cluster local position to global coords with this additional rotation added
      auto x = cluster->getLocalX() * 10;  // mm
      auto y = cluster->getLocalY() * 10;
      if (trkrId == TrkrDefs::tpcId)
      {
        y = convertTimeToZ(cluster_key, cluster);
      }

      Eigen::Vector3d clusterLocalPosition(x, y, 0);                                      // follows the convention for the acts transform of local = (x,z,y)
      Eigen::Vector3d finalCoords = overallTransformation * clusterLocalPosition / 10.0;  // convert mm back to cm

      // have to add corrections for TPC clusters after transformation to global
      if (trkrId == TrkrDefs::tpcId)
      {
        makeTpcGlobalCorrections(cluster_key, crossing, global);
      }

      if (ip == 0)
      {
        keeper(0) = (finalCoords(0) - global(0));
        keeper(1) = (finalCoords(1) - global(1));
        keeper(2) = (finalCoords(2) - global(2));
      }
      else
      {
        keeper(0) -= (finalCoords(0) - global(0));
        keeper(1) -= (finalCoords(1) - global(1));
        keeper(2) -= (finalCoords(2) - global(2));
      }

      if (m_verbosity > 1)
      {
        std::cout << "        finalCoords(0) " << finalCoords(0) << " global(0) " << global(0) << " finalCoords(1) " << finalCoords(1)
                  << " global(1) " << global(1) << " finalCoords(2) " << finalCoords(2) << " global(2) " << global(2) << std::endl;
        std::cout << "        keeper now:  keeper(0) " << keeper(0) << " keeper(1) " << keeper(1) << " keeper(2) " << keeper(2) << std::endl;
      }
    }

    // derivs vector contains:
    //   (dx/dalpha,     dy/dalpha,     dz/dalpha)     (for iangle = 0)
    //   (dx/dbeta,        dy/dbeta,        dz/dbeta)        (for iangle = 1)
    //   (dx/dgamma, dy/dgamma, dz/dgamma) (for iangle = 2)

    // Average the changes to get the estimate of the derivative
    // Check what the sign of this should be !!!!
    derivs(0) = keeper(0) / (2.0 * fabs(theseAngles[iangle]));
    derivs(1) = keeper(1) / (2.0 * fabs(theseAngles[iangle]));
    derivs(2) = keeper(2) / (2.0 * fabs(theseAngles[iangle]));
    derivs_vector.push_back(derivs);

    if (m_verbosity > 1)
    {
      std::cout << "        derivs(0) " << derivs(0) << " derivs(1) " << derivs(1) << " derivs(2) " << derivs(2) << std::endl;
    }
  }

  return derivs_vector;
}

Acts::Transform3 ActsAlignmentStates::makePerturbationTransformation(Acts::Vector3 angles)
{
  // Note: Here beta is apllied to the z axis and gamma is applied to the y axis because the geocontext transform
  // will flip those axes when transforming to global coordinates
  Eigen::AngleAxisd alpha(angles(0), Eigen::Vector3d::UnitX());
  Eigen::AngleAxisd beta(angles(2), Eigen::Vector3d::UnitY());
  Eigen::AngleAxisd gamma(angles(1), Eigen::Vector3d::UnitZ());
  Eigen::Quaternion<double> q = gamma * beta * alpha;
  Eigen::Matrix3d perturbationRotation = q.matrix();

  // combine rotation and translation into an affine matrix
  Eigen::Vector3d nullTranslation(0, 0, 0);
  Acts::Transform3 perturbationTransformation;
  perturbationTransformation.linear() = perturbationRotation;
  perturbationTransformation.translation() = nullTranslation;

  return perturbationTransformation;
}
float ActsAlignmentStates::convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster* cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = m_tGeometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89;                // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;  // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if (side == 0) zloc = -zloc;
  float z = zloc * 10.0;

  return z;
}
