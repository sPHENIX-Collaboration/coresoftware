#include "PHActsGSF.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGsfTrackFittingAlgorithm.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/BetheHeitlerApprox.hpp>

#include <TDatabasePDG.h>

//____________________________________________________________________________..
PHActsGSF::PHActsGSF(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHActsGSF::~PHActsGSF()
{
}

//____________________________________________________________________________..
int PHActsGSF::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHActsGSF::InitRun begin" << std::endl;
  }

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto bha = Acts::Experimental::makeDefaultBetheHeitlerApprox();
  ActsGsfTrackFittingAlgorithm gsf;
  m_fitCfg.fit = gsf.makeGsfFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField,
      bha, 
      4, Acts::FinalReductionMethod::eMean, true, false);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHActsGSF::process_event(PHCompositeNode*)
{
  auto logLevel = Acts::Logging::FATAL;
  if (Verbosity() > 4)
  {
    logLevel = Acts::Logging::VERBOSE;
  }

  auto logger = Acts::getDefaultLogger("PHActsGSF", logLevel);

  for (const auto& [key, track] : *m_trackMap)
  {
    auto pSurface = makePerigee(track);
    const auto seed = makeSeed(track, pSurface);

    ActsTrackFittingAlgorithm::MeasurementContainer measurements;
    TrackSeed* tpcseed = track->get_tpc_seed();
    TrackSeed* silseed = track->get_silicon_seed();

    /// We only fit full sPHENIX tracks
    if (!silseed or !tpcseed)
    {
      continue;
    }

    auto crossing = silseed->get_crossing();
    if (crossing == SHRT_MAX)
    {
      continue;
    }

    auto sourceLinks = getSourceLinks(tpcseed, measurements, crossing);
    auto silSourceLinks = getSourceLinks(silseed, measurements, crossing);

    for (auto& siSL : silSourceLinks)
    {
      sourceLinks.push_back(siSL);
    }

    /// Acts requires a wrapped vector, so we need to replace the
    /// std::vector contents with a wrapper vector to get the memory
    /// access correct
    std::vector<std::reference_wrapper<const SourceLink>> wrappedSls;
    for (const auto& sl : sourceLinks)
    {
      wrappedSls.push_back(std::cref(sl));
    }

    Calibrator calibrator(measurements);
    auto magcontext = m_tGeometry->geometry().magFieldContext;
    auto calcontext = m_tGeometry->geometry().calibContext;

    auto ppoptions = Acts::PropagatorPlainOptions();
    ppoptions.absPdgCode = m_pHypothesis;
    ppoptions.mass = TDatabasePDG::Instance()->GetParticle(
                                                 m_pHypothesis)
                         ->Mass() *
                     Acts::UnitConstants::GeV;
    ActsTrackFittingAlgorithm::GeneralFitterOptions options{
        m_tGeometry->geometry().getGeoContext(),
        magcontext,
        calcontext,
        calibrator, &(*pSurface), Acts::LoggerWrapper(*logger),
        ppoptions};
    if (Verbosity() > 2)
    {
      std::cout << "calling gsf with position "
                << seed.position(m_tGeometry->geometry().getGeoContext()).transpose()
                << " and momentum " << seed.momentum().transpose()
                << std::endl;
    }
    auto result = fitTrack(wrappedSls, seed, options);
    std::cout << "result returned" << std::endl;
    if (result.ok())
    {
      const FitResult& output = result.value();
      updateTrack(output, track);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::shared_ptr<Acts::PerigeeSurface> PHActsGSF::makePerigee(SvtxTrack* track) const
{
  SvtxVertex* vertex = m_vertexMap->get(track->get_vertex_id());

  Acts::Vector3 vertexpos(vertex->get_x() * Acts::UnitConstants::cm,
                          vertex->get_y() * Acts::UnitConstants::cm,
                          vertex->get_z() * Acts::UnitConstants::cm);

  return Acts::Surface::makeShared<Acts::PerigeeSurface>(
      vertexpos);
}

ActsTrackFittingAlgorithm::TrackParameters PHActsGSF::makeSeed(SvtxTrack* track,
                                                  std::shared_ptr<Acts::PerigeeSurface> psurf) const
{
  Acts::Vector4 fourpos(track->get_x() * Acts::UnitConstants::cm,
                        track->get_y() * Acts::UnitConstants::cm,
                        track->get_z() * Acts::UnitConstants::cm,
                        10 * Acts::UnitConstants::ns);

  int charge = track->get_charge();
  Acts::Vector3 momentum(track->get_px(),
                         track->get_py(),
                         track->get_pz());

  ActsTransformations transformer;
  auto cov = transformer.rotateSvtxTrackCovToActs(track);

  return ActsTrackFittingAlgorithm::TrackParameters::create(psurf,
                                               m_tGeometry->geometry().getGeoContext(),
                                               fourpos,
                                               momentum,
                                               charge / momentum.norm(),
                                               cov)
      .value();
}

SourceLinkVec PHActsGSF::getSourceLinks(TrackSeed* track,
                                        ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
                                        const short int& crossing)
{
  SourceLinkVec sls;
  // loop over all clusters
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;

  for (auto clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
  {
    auto key = *clusIter;
    auto cluster = m_clusterContainer->findCluster(key);
    if (!cluster)
    {
      if (Verbosity() > 0) std::cout << "Failed to get cluster with key " << key << std::endl;
      continue;
    }

    auto subsurfkey = cluster->getSubSurfKey();

    /// Make a safety check for clusters that couldn't be attached
    /// to a surface
    auto surf = m_tGeometry->maps().getSurface(key, cluster);
    if (!surf)
    {
      continue;
    }

    unsigned int trkrid = TrkrDefs::getTrkrId(key);
    unsigned int side = TpcDefs::getSide(key);

    // For the TPC, cluster z has to be corrected for the crossing z offset, distortion, and TOF z offset
    // we do this locally here and do not modify the cluster, since the cluster may be associated with multiple silicon tracks

    auto global = m_tGeometry->getGlobalPosition(key, cluster);

    if (Verbosity() > 0)
    {
      std::cout << " zinit " << global[2] << " xinit " << global[0] << " yinit " << global[1] << " side " << side << " crossing " << crossing
                << " cluskey " << key << " subsurfkey " << subsurfkey << std::endl;
    }

    if (trkrid == TrkrDefs::tpcId)
    {
      // make all corrections to global position of TPC cluster
      float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
      global[2] = z;

      // apply distortion corrections
      if (m_dccStatic)
      {
        global = m_distortionCorrection.get_corrected_position(global, m_dccStatic);
      }
      if (m_dccAverage)
      {
        global = m_distortionCorrection.get_corrected_position(global, m_dccAverage);
      }
      if (m_dccFluctuation)
      {
        global = m_distortionCorrection.get_corrected_position(global, m_dccFluctuation);
      }
    }

    // add the global positions to a vector to give to the cluster mover
    global_raw.push_back(std::make_pair(key, global));

  }  // end loop over clusters here

  // move the cluster positions back to the original readout surface
  auto global_moved = m_clusterMover.processTrack(global_raw);

  // loop over global positions returned by cluster mover
  for (int i = 0; i < global_moved.size(); ++i)
  {
    TrkrDefs::cluskey cluskey = global_moved[i].first;
    Acts::Vector3 global = global_moved[i].second;
    auto cluster = m_clusterContainer->findCluster(cluskey);

    Surface surf = m_tGeometry->maps().getSurface(cluskey, cluster);
    TrkrDefs::subsurfkey subsurfkey;

    unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);
    if (trkrid == TrkrDefs::tpcId)
    {
      // get the new surface corresponding to this global position
      TrkrDefs::hitsetkey tpcHitSetKey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
      surf = m_tGeometry->get_tpc_surface_from_coords(tpcHitSetKey,
                                                      global,
                                                      subsurfkey);
    }

    if (!surf)
    {
      continue;
    }

    // get local coordinates
    Acts::Vector2 localPos;
    Acts::Vector3 normal = surf->normal(m_tGeometry->geometry().getGeoContext());
    auto local = surf->globalToLocal(m_tGeometry->geometry().getGeoContext(),
                                     global * Acts::UnitConstants::cm,
                                     normal);

    if (local.ok())
    {
      localPos = local.value() / Acts::UnitConstants::cm;
    }
    else
    {
      /// otherwise take the manual calculation
      Acts::Vector3 loct = surf->transform(m_tGeometry->geometry().getGeoContext()).inverse() * global;
      loct /= Acts::UnitConstants::cm;

      localPos(0) = loct(0);
      localPos(1) = loct(1);
    }

    if (Verbosity() > 0)
    {
      std::cout << " cluster global after mover: " << global << std::endl;
      std::cout << " cluster local X " << cluster->getLocalX() << " cluster local Y " << cluster->getLocalY() << std::endl;
      std::cout << " new      local X " << localPos(0) << " new       local Y " << localPos(1) << std::endl;
    }

    Acts::ActsVector<2> loc;
    loc[Acts::eBoundLoc0] = localPos(0) * Acts::UnitConstants::cm;
    loc[Acts::eBoundLoc1] = localPos(1) * Acts::UnitConstants::cm;
    std::array<Acts::BoundIndices, 2> indices;
    indices[0] = Acts::BoundIndices::eBoundLoc0;
    indices[1] = Acts::BoundIndices::eBoundLoc1;
    Acts::ActsSymMatrix<2> cov = Acts::ActsSymMatrix<2>::Zero();
    if (m_cluster_version == 3)
    {
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
          cluster->getActsLocalError(0, 0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc0, Acts::eBoundLoc1) =
          cluster->getActsLocalError(0, 1) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc0) =
          cluster->getActsLocalError(1, 0) * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
          cluster->getActsLocalError(1, 1) * Acts::UnitConstants::cm2;
    }
    else if (m_cluster_version == 4)
    {
      double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
      auto para_errors = _ClusErrPara.get_cluster_error(track, cluster, clusRadius, cluskey);
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = para_errors.first * Acts::UnitConstants::cm2;
      cov(Acts::eBoundLoc0, Acts::eBoundLoc1) = 0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 0;
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = para_errors.second * Acts::UnitConstants::cm2;
    }

    ActsSourceLink::Index index = measurements.size();

    SourceLink sl(surf->geometryId(), index, cluskey);

    Acts::Measurement<Acts::BoundIndices, 2> meas(sl, indices, loc, cov);
    if (Verbosity() > 3)
    {
      std::cout << "source link " << sl.index() << ", loc : "
                << loc.transpose() << std::endl
                << ", cov : " << cov.transpose() << std::endl
                << " geo id " << sl.geometryId() << std::endl;
      std::cout << "Surface : " << std::endl;
      surf.get()->toStream(m_tGeometry->geometry().getGeoContext(), std::cout);
      std::cout << std::endl;
      std::cout << "Cluster error " << sqrt(cov(Acts::eBoundLoc0, Acts::eBoundLoc0)) << " , " << sqrt(cov(Acts::eBoundLoc1, Acts::eBoundLoc1)) << std::endl;
      std::cout << "For key " << cluskey << " with local pos " << std::endl
                << localPos(0) << ", " << localPos(1)
                << std::endl;
    }

    sls.push_back(sl);
    measurements.push_back(meas);
  }

  return sls;
}

ActsTrackFittingAlgorithm::TrackFitterResult PHActsGSF::fitTrack(
    const std::vector<std::reference_wrapper<const SourceLink>>& sourceLinks,
    const ActsTrackFittingAlgorithm::TrackParameters& seed,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& options)
{
  auto mtj = std::make_shared<Acts::VectorMultiTrajectory>();
  return (*m_fitCfg.fit)(sourceLinks, seed, options,mtj);
}

void PHActsGSF::updateTrack(const FitResult& result, SvtxTrack* track)
{
  std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
  trackTips.reserve(1);
  trackTips.emplace_back(result.lastMeasurementIndex);
  ActsExamples::Trajectories::IndexedParameters indexedParams;
  
  if (result.fittedParameters)
  {
    indexedParams.emplace(result.lastMeasurementIndex,
                          result.fittedParameters.value());
    Trajectory traj(result.fittedStates, trackTips, indexedParams);

    updateSvtxTrack(traj, track);
  }
}

void PHActsGSF::updateSvtxTrack(const Trajectory& traj, SvtxTrack* track)
{
  std::cout << "updating svtxtrack" << std::endl;
  const auto& mj = traj.multiTrajectory();
  const auto& tips = traj.tips();
  const auto& tracktip = tips.front();

  const auto& params = traj.trackParameters(tracktip);
  const auto trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(mj, tracktip);

  if (Verbosity() > 1)
  {
    std::cout << "Old track parameters: " << std::endl
              << "   (" << track->get_x()
              << ", " << track->get_y() << ", " << track->get_z()
              << ")" << std::endl
              << "   (" << track->get_px() << ", " << track->get_py()
              << ", " << track->get_pz() << ")" << std::endl;
    std::cout << "New GSF track parameters: " << std::endl
              << "   " << params.position(m_tGeometry->geometry().getGeoContext()).transpose()
              << std::endl
              << "   " << params.momentum().transpose()
              << std::endl;
  }

  /// Will create new states
  track->clear_states();

  // create a state at pathlength = 0.0
  // This state holds the track parameters, which will be updated below
  float pathlength = 0.0;
  SvtxTrackState_v1 out(pathlength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  track->insert_state(&out);

  track->set_x(params.position(m_tGeometry->geometry().getGeoContext())(0) / Acts::UnitConstants::cm);
  track->set_y(params.position(m_tGeometry->geometry().getGeoContext())(1) / Acts::UnitConstants::cm);
  track->set_z(params.position(m_tGeometry->geometry().getGeoContext())(2) / Acts::UnitConstants::cm);

  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));
  track->set_charge(params.charge());
  track->set_chisq(trajState.chi2Sum);
  track->set_ndf(trajState.NDF);

  ActsTransformations transformer;
  transformer.setVerbosity(Verbosity());

  if (params.covariance())
  {
    auto rotatedCov = transformer.rotateActsCovToSvtxTrack(params);
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        track->set_error(i, j, rotatedCov(i, j));
      }
    }
  }

  transformer.fillSvtxTrackStates(mj, tracktip, track, m_tGeometry->geometry().getGeoContext());
}

//____________________________________________________________________________..
int PHActsGSF::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHActsGSF::getNodes(PHCompositeNode* topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << PHWHERE << " The input track map is not available. Exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE << "The input cluster container is not available. Exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "The input Acts tracking geometry is not available. Exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc distortion corrections
  m_dccStatic = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dccStatic)
  {
    std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl;
  }

  m_dccAverage = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (m_dccAverage)
  {
    std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl;
  }

  m_dccFluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerFluctuation");
  if (m_dccFluctuation)
  {
    std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "Vertex map unavailable, exiting PHActsGSF" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
