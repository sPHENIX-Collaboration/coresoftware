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
  template <class T>
  inline constexpr T get_r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }
}  // namespace

void ActsAlignmentStates::fillAlignmentStateMap(const ActsTrackFittingAlgorithm::TrackContainer& tracks,
                                                const std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                                                SvtxTrack* track,
                                                const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  const auto& mj = tracks.trackStateContainer();
  const auto& trackTip = tips.front();
  const auto crossing = track->get_silicon_seed()->get_crossing();

  ActsPropagator propagator(m_tGeometry);
  if (m_fieldMap.find(".root") == std::string::npos)
  {
    propagator.constField();
    propagator.setConstFieldValue(std::stod(m_fieldMap));
  }

  SvtxAlignmentStateMap::StateVec statevec;
  if (m_verbosity > 2)
  {
    std::cout << "Beginning alignment state creation for track "
              << track->get_id() << std::endl;
  }

  std::vector<Acts::BoundIndices> indices{Acts::eBoundLoc0, Acts::eBoundLoc1, Acts::eBoundPhi, Acts::eBoundTheta, Acts::eBoundQOverP, Acts::eBoundTime};

  auto silseed = track->get_silicon_seed();
  int nmaps = 0;
  int nintt = 0;
  for (auto iter = silseed->begin_cluster_keys();
       iter != silseed->end_cluster_keys();
       ++iter)
  {
    TrkrDefs::cluskey ckey = *iter;
    if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId)
    {
      nmaps++;
    }
    if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::inttId)
    {
      nintt++;
    }
  }

  if (nmaps < 2 or nintt < 2)
  {
    return;
  }

  //! make sure the track was fully fit through the mvtx
  SvtxTrackState* firststate = (*std::next(track->begin_states(), 1)).second;
  if (get_r(firststate->get_x(), firststate->get_y()) > 5.)
  {
    return;
  }

  mj.visitBackwards(trackTip, [&](const auto& state)
                    {
    /// Collect only track states which were used in smoothing of KF and are measurements
    if (not state.hasSmoothed() or
        not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag))
    {
      return true;
    }
      
    const auto& surface = state.referenceSurface();
    auto sl = state.getUncalibratedSourceLink().template get<ActsSourceLink>();
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

    clus_sigma(2) = clus->getZError() * Acts::UnitConstants::cm;
    clus_sigma(0) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
    clus_sigma(1) = clus->getRPhiError() / sqrt(2) * Acts::UnitConstants::cm;
           
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

    Acts::Vector3 tangent = Acts::makeDirectionFromPhiTheta(phi,theta);
    //! opposite convention for coordinates in Acts
    tangent *= -1;

    if(m_verbosity > 2)
      {
	std::cout << "tangent vector to track state is " << tangent.transpose() << std::endl;
      }

    std::pair<Acts::Vector3, Acts::Vector3> projxy = 
      get_projectionXY(surface, tangent);

    Acts::Vector3 sensorCenter = surface.center(m_tGeometry->geometry().getGeoContext());
    Acts::Vector3 OM = stateGlobal - sensorCenter;

    auto globDeriv = makeGlobalDerivatives(OM, projxy);

    if(m_verbosity > 2)
      {
	std::cout << "   global deriv calcs" << std::endl
		  << "stateGlobal: " << stateGlobal.transpose()
		  << ", sensor center " << sensorCenter.transpose() << std::endl
		  << ", OM " << OM.transpose() << std::endl << "   projxy "
		  << projxy.first.transpose() << ", " 
		  << projxy.second.transpose() << std::endl
		  << "global derivatives " << std::endl << globDeriv << std::endl;
      }

    //! this is the derivative of the state wrt to Acts track parameters
    //! e.g. (d_0, z_0, phi, theta, q/p, t)
    auto localDeriv = state.effectiveProjector() * state.jacobian();
    if(m_verbosity > 2)
      {
	std::cout << "local deriv " << std::endl << localDeriv << std::endl;
      }  

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
                                           const std::pair<Acts::Vector3, Acts::Vector3>& projxy)
{
  SvtxAlignmentState::GlobalMatrix globalder = SvtxAlignmentState::GlobalMatrix::Zero();
  Acts::SquareMatrix3 unit = Acts::SquareMatrix3::Identity();

  //! x residual rotations
  globalder(0, 0) = ((unit.col(0)).cross(OM)).dot(projxy.first);
  globalder(0, 1) = ((unit.col(1)).cross(OM)).dot(projxy.first);
  globalder(0, 2) = ((unit.col(2)).cross(OM)).dot(projxy.first);
  //! x residual translations
  globalder(0, 3) = unit.col(0).dot(projxy.first);
  globalder(0, 4) = unit.col(1).dot(projxy.first);
  globalder(0, 5) = unit.col(2).dot(projxy.first);

  //! y residual rotations
  globalder(1, 0) = ((unit.col(0)).cross(OM)).dot(projxy.second);
  globalder(1, 1) = ((unit.col(1)).cross(OM)).dot(projxy.second);
  globalder(1, 2) = ((unit.col(2)).cross(OM)).dot(projxy.second);
  //! y residual translations
  globalder(1, 3) = unit.col(0).dot(projxy.second);
  globalder(1, 4) = unit.col(1).dot(projxy.second);
  globalder(1, 5) = unit.col(2).dot(projxy.second);

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
  Acts::Vector3 xloc(1.0, 0.0, 0.0);
  Acts::Vector3 xglob = surface.transform(m_tGeometry->geometry().getGeoContext()) * xloc;

  Acts::Vector3 yloc(0.0, 1.0, 0.0);
  Acts::Vector3 yglob = surface.transform(m_tGeometry->geometry().getGeoContext()) * yloc;

  Acts::Vector3 X = (xglob - sensorCenter) / (xglob - sensorCenter).norm();
  Acts::Vector3 Y = (yglob - sensorCenter) / (yglob - sensorCenter).norm();

  projx = X - (tangent.dot(X) / tangent.dot(Z)) * Z;
  projy = Y - (tangent.dot(Y) / tangent.dot(Z)) * Z;
  if (m_verbosity > 2)
    std::cout << "projxy " << projx.transpose() << ", " << projy.transpose() << std::endl;
  return std::make_pair(projx, projy);
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
