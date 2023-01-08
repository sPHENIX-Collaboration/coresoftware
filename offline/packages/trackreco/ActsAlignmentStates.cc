#include "ActsAlignmentStates.h"

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

void ActsAlignmentStates::fillAlignmentStateMap(Trajectory traj,
                                                SvtxTrack* track)
{
  const auto mj = traj.multiTrajectory();
  const auto& tips = traj.tips();
  const auto& trackTip = tips.front();
  const auto crossing = track->get_silicon_seed()->get_crossing();
  SvtxAlignmentStateMap::StateVec statevec;

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

    if (m_verbosity > 2)
      {
	std::cout << "sl index and ckey " << sl.index() << ", "
		  << sl.cluskey() << std::endl;
      }

    auto clus = m_clusterMap->findCluster(ckey);
    const auto trkrId = TrkrDefs::getTrkrId(ckey);

    /// Gets the global parameters from the state
    const Acts::FreeVector freeParams =
        Acts::MultiTrajectoryHelpers::freeSmoothed(m_tGeometry->geometry().getGeoContext(), state);
 
    /// Calculate the residual in global coordinates
    Acts::Vector3 clusGlobal = m_tGeometry->getGlobalPosition(ckey, clus);
    if (trkrId == TrkrDefs::tpcId)
    {
      makeTpcGlobalCorrections(ckey, crossing, clusGlobal);
    }

    /// convert to acts units
    clusGlobal *= Acts::UnitConstants::cm;

    const Acts::FreeVector globalStateParams = Acts::detail::transformBoundToFreeParameters(surface, m_tGeometry->geometry().getGeoContext(), state.smoothed());
    Acts::Vector3 stateGlobal = globalStateParams.segment<3>(Acts::eFreePos0);
    Acts::Vector3 residual = clusGlobal - stateGlobal;

    Acts::Vector3 clus_sigma(0, 0, 0);
      if (m_clusterVersion == 3)
      {
        clus_sigma(2) = clus->getZError() * Acts::UnitConstants::cm;
        clus_sigma(0) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
        clus_sigma(1) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
      }
      else if (m_clusterVersion == 4)
      {
        double clusRadius = sqrt(clusGlobal(0) * clusGlobal(0) + clusGlobal(1) * clusGlobal(1));
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
                << "state global is " << stateGlobal.transpose() << std::endl
                << "Residual is " << residual.transpose() << std::endl;
      std::cout << "   clus errors are " << clus_sigma.transpose() << std::endl;
    }

    // Get the derivative of alignment (global) parameters w.r.t. measurement or residual
    const Acts::Vector3 direction = freeParams.segment<3>(Acts::eFreeDir0);
    // The derivative of free parameters w.r.t. path length. @note Here, we
    // assumes a linear track model, i.e. negecting the change of track
    // direction. Otherwise, we need to know the magnetic field at the free
    // parameters
    Acts::FreeVector pathDerivative = Acts::FreeVector::Zero();
    pathDerivative.head<3>() = direction;

    const Acts::ActsDynamicMatrix H = state.effectiveProjector();

    /// Acts residual, in local coordinates
    //auto actslocres = state.effectiveCalibrated() - H * state.smoothed();

    // Get the derivative of bound parameters w.r.t. alignment parameters
    Acts::AlignmentToBoundMatrix d =
        surface.alignmentToBoundDerivative(m_tGeometry->geometry().getGeoContext(), freeParams, pathDerivative);
    // Get the derivative of bound parameters wrt track parameters
    Acts::FreeToBoundMatrix j = surface.freeToBoundJacobian(m_tGeometry->geometry().getGeoContext(), freeParams);

    // derivative of residual wrt track parameters
    auto dLocResTrack = -H * j;
    // derivative of residual wrt alignment parameters
    auto dLocResAlignment = -H * d;

    if (m_verbosity > 4)
    {
      std::cout << " derivative of resiudal wrt track params " << std::endl
                << dLocResTrack << std::endl
                << " derivative of residual wrt alignment params " << std::endl
                << dLocResAlignment << std::endl;
    }

    /// The above matrices are in Acts local coordinates, so we need to convert them to
    /// global coordinates with a rotation d(l0,l1)/d(x,y,z)

    const double cosTheta = direction.z();
    const double sinTheta = std::sqrt(square(direction.x()) + square(direction.y()));
    const double invSinTheta = 1. / sinTheta;
    const double cosPhi = direction.x() * invSinTheta;
    const double sinPhi = direction.y() * invSinTheta;
    Acts::ActsMatrix<3, 2> rot = Acts::ActsMatrix<3, 2>::Zero();
    rot(0, 0) = -sinPhi;
    rot(0, 1) = -cosPhi * cosTheta;
    rot(1, 0) = cosPhi;
    rot(1, 1) = -sinPhi * cosTheta;
    rot(2, 1) = sinTheta;

    /// this is a 3x6 matrix now of d(x_res,y_res,z_res)/d(dx,dy,dz,alpha,beta,gamma)
    const auto dGlobResAlignment = rot * dLocResAlignment;

    /// Acts global coordinates are (x,y,z,t,hatx,haty,hatz,q/p)
    /// So now rotate to x,y,z, px,py,pz
    const double p = 1. / abs(globalStateParams[Acts::eFreeQOverP]);
    const double p2 = square(p);
    Acts::ActsMatrix<6, 8> sphenixRot = Acts::ActsMatrix<6, 8>::Zero();
    sphenixRot(0, 0) = 1;
    sphenixRot(1, 1) = 1;
    sphenixRot(2, 2) = 1;
    sphenixRot(3, 4) = p;
    sphenixRot(4, 5) = p;
    sphenixRot(5, 6) = p;
    sphenixRot(3, 7) = direction.x() * p2;
    sphenixRot(4, 7) = direction.y() * p2;
    sphenixRot(5, 7) = direction.z() * p2;

    const auto dGlobResTrack = rot * dLocResTrack * sphenixRot.transpose();

    if (m_verbosity > 3)
    {
      std::cout << "derivative of residual wrt alignment parameters glob " << std::endl
                << dGlobResAlignment << std::endl;
      std::cout << "derivative of residual wrt track parameters glob " << std::endl
                << dGlobResTrack << std::endl;
    }


    auto svtxstate = std::make_unique<SvtxAlignmentState_v1>();
    svtxstate->set_residual(residual);
    svtxstate->set_local_derivative_matrix(dGlobResTrack);
    svtxstate->set_global_derivative_matrix(dGlobResAlignment);
    svtxstate->set_cluster_key(ckey);
    
    if(m_analytic)
      {
	auto surf = m_tGeometry->maps().getSurface(ckey, clus);
	auto anglederivs = getDerivativesAlignmentAngles(clusGlobal, ckey, 
							 clus, surf, 
							 crossing);
        SvtxAlignmentState::GlobalMatrix analytic = SvtxAlignmentState::GlobalMatrix::Zero();

	getGlobalDerivatives(anglederivs,analytic);
       	svtxstate->set_global_derivative_matrix(analytic);
      }

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
