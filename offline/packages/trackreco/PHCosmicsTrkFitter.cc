#include "PHCosmicsTrkFitter.h"

/// Tracking includes
#include <trackbase/Calibrator.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentStateMap_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <micromegas/MicromegasDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include <TDatabasePDG.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  // check vector validity
  inline bool is_valid(const Acts::Vector3 vec)
  {
    return !(std::isnan(vec.x()) || std::isnan(vec.y()) || std::isnan(vec.z()));
  }
  template <class T>
  inline T square(const T& x)
  {
    return x * x;
  }
}  // namespace

#include <trackbase/alignmentTransformationContainer.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

PHCosmicsTrkFitter::PHCosmicsTrkFitter(const std::string& name)
  : SubsysReco(name)
  , m_trajectories(nullptr)
{
}

int PHCosmicsTrkFitter::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Setup PHCosmicsTrkFitter" << std::endl;
  }

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_alignStates.distortionContainers(_dcc_static, _dcc_average, _dcc_fluctuation);
  m_alignStates.actsGeometry(m_tGeometry);
  m_alignStates.clusters(m_clusterContainer);
  m_alignStates.stateMap(m_alignmentStateMap);
  m_alignStates.verbosity(Verbosity());
  m_alignStates.fieldMap(m_fieldMap);

  m_fitCfg.fit = ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField);

  m_outlierFinder.verbosity = Verbosity();
  std::map<long unsigned int, float> chi2Cuts;
  chi2Cuts.insert(std::make_pair(10, 4));
  chi2Cuts.insert(std::make_pair(12, 4));
  chi2Cuts.insert(std::make_pair(14, 9));
  chi2Cuts.insert(std::make_pair(16, 4));
  m_outlierFinder.chi2Cuts = chi2Cuts;
  if (m_useOutlierFinder)
  {
    m_fitCfg.fit->outlierFinder(m_outlierFinder);
  }

  auto cellgeo =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  if (cellgeo)
  {
    _clusterMover.initialize_geometry(cellgeo);
  }

  if (m_actsEvaluator)
  {
    m_evaluator = std::make_unique<ActsEvaluator>(m_evalname);
    m_evaluator->Init(topNode);
    m_evaluator->verbosity(Verbosity());
  }

  if (Verbosity() > 1)
  {
    std::cout << "Finish PHCosmicsTrkFitter Setup" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsTrkFitter::process_event(PHCompositeNode* topNode)
{
  m_event++;

  auto logLevel = Acts::Logging::FATAL;

  if (m_actsEvaluator)
  {
    m_evaluator->next_event(topNode);
  }

  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << "Events processed: " << m_event << std::endl;
    std::cout << "Start PHCosmicsTrkFitter::process_event" << std::endl;
    if (Verbosity() > 30)
      logLevel = Acts::Logging::VERBOSE;
  }

  /// Fill an additional track map if using the acts evaluator
  /// for proto track comparison to fitted track
  if (m_actsEvaluator)
  {
    /// wipe at the beginning of every new fit pass, so that the seeds
    /// are whatever is currently in SvtxTrackMap
    m_seedTracks->clear();
    for (const auto& [key, track] : *m_trackMap)
    {
      m_seedTracks->insert(track);
    }
  }

  loopTracks(logLevel);

  // put this in the output file
  if (Verbosity() > 0)
  {
    std::cout << " SvtxTrackMap size is now " << m_trackMap->size()
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsTrkFitter::ResetEvent(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "Reset PHCosmicsTrkFitter" << std::endl;
  }

  m_trajectories->clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsTrkFitter::End(PHCompositeNode* /*topNode*/)
{
  if (m_actsEvaluator)
  {
    m_evaluator->End();
  }

  if (Verbosity() > 0)
  {
    std::cout << "The Acts track fitter had " << m_nBadFits
              << " fits return an error" << std::endl;

    std::cout << "Finished PHCosmicsTrkFitter" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHCosmicsTrkFitter::loopTracks(Acts::Logging::Level logLevel)
{
  auto logger = Acts::getDefaultLogger("PHCosmicsTrkFitter", logLevel);

  if (Verbosity() > 0)
  {
    std::cout << " seed map size " << m_seedMap->size() << std::endl;
  }

  for (auto trackiter = m_seedMap->begin(); trackiter != m_seedMap->end();
       ++trackiter)
  {
    TrackSeed* track = *trackiter;
    if (!track)
    {
      continue;
    }

    unsigned int tpcid = track->get_tpc_seed_index();
    unsigned int siid = track->get_silicon_seed_index();

    // get the crossing number
    auto siseed = m_siliconSeeds->get(siid);
    short crossing = 0;

    auto tpcseed = m_tpcSeeds->get(tpcid);
    if (Verbosity() > 1)
    {
      std::cout << "TPC id " << tpcid << std::endl;
      std::cout << "Silicon id " << siid << std::endl;
    }

    /// Need to also check that the tpc seed wasn't removed by the ghost finder
    if (!tpcseed)
    {
      continue;
    }

    if (Verbosity() > 0)
    {
      if (siseed) std::cout << " silicon seed position is (x,y,z) = " << siseed->get_x() << "  " << siseed->get_y() << "  " << siseed->get_z() << std::endl;
      std::cout << " tpc seed position is (x,y,z) = " << tpcseed->get_x() << "  " << tpcseed->get_y() << "  " << tpcseed->get_z() << std::endl;
    }

    ActsTrackFittingAlgorithm::MeasurementContainer measurements;
    int charge = 0;
    float cosmicslope = 0;
    SourceLinkVec sourceLinks;
    if (siseed) sourceLinks = getSourceLinks(siseed, measurements, crossing, charge, cosmicslope);
    const auto tpcSourceLinks = getSourceLinks(tpcseed, measurements, crossing, charge, cosmicslope);

    sourceLinks.insert(sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end());

    tpcseed->circleFitByTaubin(m_clusterContainer, m_tGeometry, 0, 58);

    float tpcR = fabs(1. / tpcseed->get_qOverR());
    float tpcx = tpcseed->get_X0();
    float tpcy = tpcseed->get_Y0();

    const auto intersect =
        TrackFitUtils::circle_circle_intersection(m_vertexRadius,
                                                  tpcR, tpcx, tpcy);
    float intx, inty;

    if (std::get<1>(intersect) < std::get<3>(intersect))
    {
      intx = std::get<0>(intersect);
      inty = std::get<1>(intersect);
    }
    else
    {
      intx = std::get<2>(intersect);
      inty = std::get<3>(intersect);
    }

    float slope = tpcseed->get_slope();
    float intz = m_vertexRadius * slope + tpcseed->get_Z0();

    Acts::Vector3 inter(intx, inty, intz);

    std::vector<float> tpcparams{tpcR, tpcx, tpcy, tpcseed->get_slope(),
                                 tpcseed->get_Z0()};
    auto tangent = TrackFitUtils::get_helix_tangent(tpcparams,
                                                    inter);

    auto tan = tangent.second;
    auto pca = tangent.first;

    float p;
    if (m_fieldMap.find(".root") != std::string::npos)
    {
      p = tpcseed->get_p();
    }
    else
    {
      p = cosh(tpcseed->get_eta()) * fabs(1. / tpcseed->get_qOverR()) * (0.3 / 100) * std::stod(m_fieldMap);
    }

    tan *= p;

    //! if we got the opposite seed then z will be backwards, so we take the
    //! value of tan.z() multiplied by the sign of the slope determined for
    //! the full cosmic track
    //! same with px/py since a single cosmic produces two seeds that bend
    //! in opposite directions

    Acts::Vector3 momentum(charge < 0 ? tan.x() : tan.x() * -1,
                           charge < 0 ? tan.y() : tan.y() * -1,
                           cosmicslope > 0 ? fabs(tan.z()) : -1 * fabs(tan.z()));
    Acts::Vector3 position(pca.x(), pca.y(),
                           momentum.z() > 0 ? (slope < 0 ? intz : m_vertexRadius * slope * -1 + tpcseed->get_Z0()) : (slope > 0 ? intz : m_vertexRadius * slope * -1 + tpcseed->get_Z0()));

    position *= Acts::UnitConstants::cm;
    if (!is_valid(momentum)) continue;

    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        Acts::Vector3(0, -1 * m_vertexRadius * Acts::UnitConstants::cm, 0));
    auto actsFourPos = Acts::Vector4(position(0), position(1),
                                     position(2),
                                     10 * Acts::UnitConstants::ns);

    Acts::BoundSymMatrix cov = setDefaultCovariance();

    //! Acts requires a wrapped vector, so we need to replace the
    //! std::vector contents with a wrapper vector to get the memory
    //! access correct
    std::vector<Acts::SourceLink> wrappedSls;
    for (const auto& sl : sourceLinks)
    {
      wrappedSls.push_back(Acts::SourceLink{sl});
    }

    //! Reset the track seed with the dummy covariance
    auto seed = ActsTrackFittingAlgorithm::TrackParameters::create(
                    pSurface,
                    m_tGeometry->geometry().getGeoContext(),
                    actsFourPos,
                    momentum,
                    charge / momentum.norm(),
                    cov)
                    .value();

    if (Verbosity() > 2)
    {
      printTrackSeed(seed);
    }

    //! Set host of propagator options for Acts to do e.g. material integration
    Acts::PropagatorPlainOptions ppPlainOptions;
    ppPlainOptions.absPdgCode = m_pHypothesis;
    ppPlainOptions.mass = TDatabasePDG::Instance()->GetParticle(
                                                      m_pHypothesis)
                              ->Mass() *
                          Acts::UnitConstants::GeV;

    Calibrator calibrator{measurements};

    auto magcontext = m_tGeometry->geometry().magFieldContext;
    auto calibcontext = m_tGeometry->geometry().calibContext;

    ActsTrackFittingAlgorithm::GeneralFitterOptions
        kfOptions{
            m_tGeometry->geometry().getGeoContext(),
            magcontext,
            calibcontext,
            calibrator,
            &(*pSurface),
	    ppPlainOptions};

    auto trackContainer = 
	std::make_shared<Acts::VectorTrackContainer>();
      auto trackStateContainer = 
	std::make_shared<Acts::VectorMultiTrajectory>();
      ActsTrackFittingAlgorithm::TrackContainer 
	tracks(trackContainer, trackStateContainer);
    auto result = fitTrack(wrappedSls, seed, kfOptions, tracks);

    /// Check that the track fit result did not return an error
    if (result.ok())
    {
      SvtxTrack_v4 newTrack;
      newTrack.set_tpc_seed(tpcseed);
      newTrack.set_crossing(crossing);
      newTrack.set_silicon_seed(siseed);

      unsigned int trid = m_trackMap->size();
      newTrack.set_id(trid);

      if (getTrackFitResult(result, track, &newTrack, tracks, measurements))
      {
        m_trackMap->insertWithKey(&newTrack, trid);
      }
    }
    else
    {
      m_nBadFits++;
      if (Verbosity() > 1)
      {
        std::cout << "Track fit failed for track " << m_seedMap->find(track)
                  << " with Acts error message "
                  << result.error() << ", " << result.error().message()
                  << std::endl;
      }
    }
  }

  return;
}

//___________________________________________________________________________________
SourceLinkVec PHCosmicsTrkFitter::getSourceLinks(
    TrackSeed* track,
    ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
    short int crossing,
    int& charge,
    float& cosmicslope)
{
  SourceLinkVec sourcelinks;

  // loop over all clusters
  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;
  int i = 0;

  for (auto clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
  {
    i++;
    auto key = *clusIter;
    auto cluster = m_clusterContainer->findCluster(key);
    if (!cluster)
    {
      if (Verbosity() > 0)
        std::cout << "Failed to get cluster with key " << key << " for track " << m_seedMap->find(track) << std::endl;
      else
        std::cout << "PHCosmicsTrkFitter :: Key: " << key << " for track " << m_seedMap->find(track) << std::endl;
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

    const unsigned int trkrid = TrkrDefs::getTrkrId(key);
    const unsigned int side = TpcDefs::getSide(key);

    // For the TPC, cluster z has to be corrected for the crossing z offset, distortion, and TOF z offset
    // we do this locally here and do not modify the cluster, since the cluster may be associated with multiple silicon tracks
    Acts::Vector3 global = m_tGeometry->getGlobalPosition(key, cluster);

    if (trkrid == TrkrDefs::tpcId)
    {
      // make all corrections to global position of TPC cluster
      float z = m_clusterCrossingCorrection.correctZ(global[2], side, crossing);
      global[2] = z;

      // apply distortion corrections
      if (_dcc_static)
      {
        global = _distortionCorrection.get_corrected_position(global, _dcc_static);
      }
      if (_dcc_average)
      {
        global = _distortionCorrection.get_corrected_position(global, _dcc_average);
      }
      if (_dcc_fluctuation)
      {
        global = _distortionCorrection.get_corrected_position(global, _dcc_fluctuation);
      }
    }

    if (Verbosity() > 0)
    {
      std::cout << " zinit " << global[2] << " xinit " << global[0] << " yinit " << global[1] << " side " << side << " crossing " << crossing
                << " cluskey " << key << " subsurfkey " << subsurfkey << std::endl;
    }

    // add the global positions to a vector to give to the cluster mover
    global_raw.push_back(std::make_pair(key, global));

  }  // end loop over clusters here

  // move the cluster positions back to the original readout surface
  auto global_moved = _clusterMover.processTrack(global_raw);

  Acts::Vector3 globalMostOuter;
  Acts::Vector3 globalSecondMostOuter(0, 999999, 0);
  float largestR = 0;
  // loop over global positions returned by cluster mover
  for (int i = 0; i < global_moved.size(); ++i)
  {
    TrkrDefs::cluskey cluskey = global_moved[i].first;
    Acts::Vector3 global = global_moved[i].second;
    float r = std::sqrt(square(global.x()) + square(global.y()));

    /// use the bottom hemisphere to determine the charge
    if (r > largestR && global.y() < 0)
    {
      globalMostOuter = global_moved[i].second;
      largestR = r;
    }

    if (m_ignoreLayer.find(TrkrDefs::getLayer(cluskey)) != m_ignoreLayer.end())
    {
      if (Verbosity() > 3)
      {
        std::cout << PHWHERE << "skipping cluster in layer "
                  << (unsigned int) TrkrDefs::getLayer(cluskey) << std::endl;
      }
      continue;
    }

    auto cluster = m_clusterContainer->findCluster(cluskey);
    Surface surf = m_tGeometry->maps().getSurface(cluskey, cluster);

    // if this is a TPC cluster, the crossing correction may have moved it across the central membrane, check the surface
    auto trkrid = TrkrDefs::getTrkrId(cluskey);
    if (trkrid == TrkrDefs::tpcId)
    {
      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
      TrkrDefs::subsurfkey new_subsurfkey = 0;
      surf = m_tGeometry->get_tpc_surface_from_coords(hitsetkey, global, new_subsurfkey);
    }

    if (!surf)
    {
      continue;
    }

    // get local coordinates
    Acts::Vector2 localPos;
    global *= Acts::UnitConstants::cm;

    Acts::Vector3 normal = surf->normal(m_tGeometry->geometry().getGeoContext());
    auto local = surf->globalToLocal(m_tGeometry->geometry().getGeoContext(),
                                     global, normal);

    if (local.ok())
    {
      localPos = local.value() / Acts::UnitConstants::cm;
    }
    else
    {
      /// otherwise take the manual calculation for the TPC
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

    double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
    auto para_errors = _ClusErrPara.get_clusterv5_modified_error(cluster, clusRadius, cluskey);
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = para_errors.first * Acts::UnitConstants::cm2;
    cov(Acts::eBoundLoc0, Acts::eBoundLoc1) = 0;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 0;
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = para_errors.second * Acts::UnitConstants::cm2;

    ActsSourceLink::Index index = measurements.size();

    SourceLink sl(surf->geometryId(), index, cluskey);
    Acts::SourceLink actsSL{sl.geometryId(), sl};
    Acts::Measurement<Acts::BoundIndices, 2> meas(std::move(actsSL), indices, loc, cov);
    if (Verbosity() > 3)
    {
      std::cout << "source link " << sl.index() << ", loc : "
                << loc.transpose() << std::endl
                << ", cov : " << cov.transpose() << std::endl
                << " geo id " << sl.geometryId() << std::endl;
      std::cout << "Surface : " << std::endl;
      surf.get()->toStream(m_tGeometry->geometry().getGeoContext(), std::cout);
      std::cout << std::endl;
      std::cout << "Cluster error " << cluster->getRPhiError() << " , " << cluster->getZError() << std::endl;
      std::cout << "For key " << cluskey << " with local pos " << std::endl
                << localPos(0) << ", " << localPos(1)
                << std::endl;
    }

    sourcelinks.push_back(sl);
    measurements.push_back(meas);
  }

  //! find the closest cluster to the outermost cluster
  float maxdr = std::numeric_limits<float>::max();
  for (int i = 0; i < global_moved.size(); i++)
  {
    if (global_moved[i].second.y() > 0) continue;

    float dr = std::sqrt(square(globalMostOuter.x()) + square(globalMostOuter.y())) - std::sqrt(square(global_moved[i].second.x()) + square(global_moved[i].second.y()));
    //! Place a dr cut to get maximum bend due to TPC clusters having
    //! larger fluctuations
    if (dr < maxdr && dr > 10)
    {
      maxdr = dr;
      globalSecondMostOuter = global_moved[i].second;
    }
  }

  //! we have to calculate phi WRT the vertex position outside the detector,
  //! not at (0,0)
  Acts::Vector3 vertex(0, -1 * m_vertexRadius, 0);
  globalMostOuter -= vertex;
  globalSecondMostOuter -= vertex;

  const auto firstphi = atan2(globalMostOuter.y(), globalMostOuter.x());
  const auto secondphi = atan2(globalSecondMostOuter.y(),
                               globalSecondMostOuter.x());
  auto dphi = secondphi - firstphi;

  if (dphi > M_PI) dphi = 2. * M_PI - dphi;
  if (dphi < -M_PI) dphi = 2 * M_PI + dphi;

  if (dphi > 0)
  {
    charge = -1;
  }
  else
  {
    charge = 1;
  }

  float r1 = std::sqrt(square(globalMostOuter.x()) + square(globalMostOuter.y()));
  float r2 = std::sqrt(square(globalSecondMostOuter.x()) + square(globalSecondMostOuter.y()));
  float z1 = globalMostOuter.z();
  float z2 = globalSecondMostOuter.z();

  cosmicslope = (r2 - r1) / (z2 - z1);

  return sourcelinks;
}

bool PHCosmicsTrkFitter::getTrackFitResult(FitResult& fitOutput, 
					   TrackSeed* seed, SvtxTrack* track, 
					   ActsTrackFittingAlgorithm::TrackContainer& tracks,
					   const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  /// Make a trajectory state for storage, which conforms to Acts track fit
  /// analysis tool
  auto& outtrack = fitOutput.value();
  std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
  trackTips.reserve(1);
  trackTips.emplace_back(outtrack.tipIndex());
  Trajectory::IndexedParameters indexedParams;

   indexedParams.emplace(std::pair{outtrack.tipIndex(),
	 ActsExamples::TrackParameters{outtrack.referenceSurface().getSharedPtr(),
	   outtrack.parameters(), outtrack.covariance()}});
  
 if (Verbosity() > 2)
    {
      std::cout << "Fitted parameters for track" << std::endl;
      std::cout << " position : " << outtrack.referenceSurface().localToGlobal(m_tGeometry->geometry().getGeoContext(), Acts::Vector2(outtrack.loc0(), outtrack.loc1()), Acts::Vector3(1,1,1)).transpose()
	
		<< std::endl;
      std::cout << "charge: "<<outtrack.charge()<<std::endl;
      std::cout << " momentum : " << outtrack.momentum().transpose()
		<< std::endl;
      std::cout << "For trackTip == " << outtrack.tipIndex() << std::endl;
    }  
  
    
  Trajectory trajectory(tracks.trackStateContainer(),
			trackTips, indexedParams);

  m_trajectories->insert(std::make_pair(track->get_id(), trajectory));

  /// Get position, momentum from the Acts output. Update the values of
  /// the proto track
  updateSvtxTrack(trackTips, indexedParams, tracks, track);
  

  if (m_commissioning)
  {
    if (track->get_silicon_seed() && track->get_tpc_seed())
    {
      m_alignStates.fillAlignmentStateMap(tracks, trackTips,
					  track, measurements);
    }
  }

  if (m_actsEvaluator)
  {
    m_evaluator->evaluateTrackFit(tracks,trackTips,indexedParams, track,
                                  seed, measurements);
  }

  return true;
}

inline ActsTrackFittingAlgorithm::TrackFitterResult PHCosmicsTrkFitter::fitTrack(
    const std::vector<Acts::SourceLink>& sourceLinks,
    const ActsTrackFittingAlgorithm::TrackParameters& seed,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& kfOptions,
    ActsTrackFittingAlgorithm::TrackContainer& tracks)
{
  return (*m_fitCfg.fit)(sourceLinks, seed, kfOptions, tracks);
}

void PHCosmicsTrkFitter::updateSvtxTrack(
     std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
     Trajectory::IndexedParameters& paramsMap,
     ActsTrackFittingAlgorithm::TrackContainer& tracks, 
     SvtxTrack* track)
{
  const auto& mj = tracks.trackStateContainer();

  /// only one track tip in the track fit Trajectory
  auto& trackTip = tips.front();

  if (Verbosity() > 2)
  {
    std::cout << "Identify (proto) track before updating with acts results " << std::endl;
    track->identify();
  }

  // create a state at pathlength = 0.0
  // This state holds the track parameters, which will be updated below
  float pathlength = 0.0;
  SvtxTrackState_v1 out(pathlength);
  out.set_x(0.0);
  out.set_y(0.0);
  out.set_z(0.0);
  track->insert_state(&out);

  auto trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

  const auto& params = paramsMap.find(trackTip)->second;

  /// Acts default unit is mm. So convert to cm
  track->set_x(params.position(m_tGeometry->geometry().getGeoContext())(0) / Acts::UnitConstants::cm);
  track->set_y(params.position(m_tGeometry->geometry().getGeoContext())(1) / Acts::UnitConstants::cm);
  track->set_z(params.position(m_tGeometry->geometry().getGeoContext())(2) / Acts::UnitConstants::cm);

  track->set_px(params.momentum()(0));
  track->set_py(params.momentum()(1));
  track->set_pz(params.momentum()(2));

  track->set_charge(params.charge());
  track->set_chisq(trajState.chi2Sum);
  track->set_ndf(trajState.NDF);

  ActsTransformations rotater;
  rotater.setVerbosity(Verbosity());

  if (params.covariance())
  {
    auto rotatedCov = rotater.rotateActsCovToSvtxTrack(params);

    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        track->set_error(i, j, rotatedCov(i, j));
      }
    }
  }

  // Also need to update the state list and cluster ID list for all measurements associated with the acts track
  // loop over acts track states, copy over to SvtxTrackStates, and add to SvtxTrack

  if (m_fillSvtxTrackStates)
  {
    rotater.fillSvtxTrackStates(mj, trackTip, track,
                                m_tGeometry->geometry().getGeoContext());
  }

  if (Verbosity() > 2)
  {
    std::cout << " Identify fitted track after updating track states:"
              << std::endl;
    track->identify();
  }

  return;
}

Acts::BoundSymMatrix PHCosmicsTrkFitter::setDefaultCovariance() const
{
  Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();

  /// Acts cares about the track covariance as it helps the KF
  /// know whether or not to trust the initial track seed or not.
  /// We reset it here to some loose values as it helps Acts improve
  /// the fitting.
  /// If the covariance is too loose, it won't be able to propagate,
  /// but if it is too tight, it will just "believe" the track seed over
  /// the hit data

  double sigmaD0 = 300 * Acts::UnitConstants::um;
  double sigmaZ0 = 300 * Acts::UnitConstants::um;
  double sigmaPhi = 1 * Acts::UnitConstants::degree;
  double sigmaTheta = 1 * Acts::UnitConstants::degree;
  double sigmaT = 1. * Acts::UnitConstants::ns;

  cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = sigmaD0 * sigmaD0;
  cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = sigmaZ0 * sigmaZ0;
  cov(Acts::eBoundTime, Acts::eBoundTime) = sigmaT * sigmaT;
  cov(Acts::eBoundPhi, Acts::eBoundPhi) = sigmaPhi * sigmaPhi;
  cov(Acts::eBoundTheta, Acts::eBoundTheta) = sigmaTheta * sigmaTheta;
  /// Acts takes this value very seriously - tuned to be in a "sweet spot"
  cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = 0.0001;

  return cov;
}

void PHCosmicsTrkFitter::printTrackSeed(const ActsTrackFittingAlgorithm::TrackParameters& seed) const
{
  std::cout
      << PHWHERE
      << " Processing proto track:"
      << std::endl;

  std::cout
      << "position: " << seed.position(m_tGeometry->geometry().getGeoContext()).transpose()
      << std::endl
      << "momentum: " << seed.momentum().transpose()
      << std::endl;

  std::cout << "charge : " << seed.charge() << std::endl;
  std::cout << "absolutemom : " << seed.absoluteMomentum() << std::endl;
}

int PHCosmicsTrkFitter::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHCosmicsTrkFitter::createNodes");
  }

  PHNodeIterator dstIter(topNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_trajectories = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsTrajectories");
  if (!m_trajectories)
  {
    m_trajectories = new std::map<const unsigned int, Trajectory>;
    auto node =
        new PHDataNode<std::map<const unsigned int, Trajectory>>(m_trajectories, "ActsTrajectories");
    svtxNode->addNode(node);
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);

  if (!m_trackMap)
  {
    m_trackMap = new SvtxTrackMap_v2;
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_trackMap, _track_map_name, "PHObject");
    svtxNode->addNode(node);
  }

  m_alignmentStateMap = findNode::getClass<SvtxAlignmentStateMap>(topNode, "SvtxAlignmentStateMap");
  if (!m_alignmentStateMap)
  {
    m_alignmentStateMap = new SvtxAlignmentStateMap_v1;
    auto node = new PHDataNode<SvtxAlignmentStateMap>(m_alignmentStateMap, "SvtxAlignmentStateMap", "PHObject");
    svtxNode->addNode(node);
  }

  if (m_actsEvaluator)
  {
    m_seedTracks = findNode::getClass<SvtxTrackMap>(topNode, _seed_track_map_name);

    if (!m_seedTracks)
    {
      m_seedTracks = new SvtxTrackMap_v2;

      PHIODataNode<PHObject>* seedNode =
          new PHIODataNode<PHObject>(m_seedTracks, _seed_track_map_name, "PHObject");
      svtxNode->addNode(seedNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsTrkFitter::getNodes(PHCompositeNode* topNode)
{
  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!m_tpcSeeds)
  {
    std::cout << PHWHERE << "TpcTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!m_siliconSeeds)
  {
    std::cout << PHWHERE << "SiliconTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE
              << "No trkr cluster container, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_seedMap = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if (!m_seedMap)
  {
    std::cout << "No Svtx seed map on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc distortion corrections
  _dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (_dcc_static)
  {
    std::cout << PHWHERE << "  found static TPC distortion correction container" << std::endl;
  }
  _dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerAverage");
  if (_dcc_average)
  {
    std::cout << PHWHERE << "  found average TPC distortion correction container" << std::endl;
  }
  _dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerFluctuation");
  if (_dcc_fluctuation)
  {
    std::cout << PHWHERE << "  found fluctuation TPC distortion correction container" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
