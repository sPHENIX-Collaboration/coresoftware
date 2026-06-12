#include "PHCosmicsTrkFitter.h"
#include "MakeSourceLinks.h"

/// Tracking includes
#include <trackbase/Calibrator.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentStateMap_v1.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v1.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <g4detectors/PHG4TpcGeomContainer.h>

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

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  // check vector validity
  inline bool is_valid(const Acts::Vector3& vec)
  {
    return !(std::isnan(vec.x()) || std::isnan(vec.y()) || std::isnan(vec.z()));
  }
  template <class T>
  inline T square(const T& x)
  {
    return x * x;
  }

  /// get radius from coordinates
  template <class T>
  T radius(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

}  // namespace

#include <trackbase/alignmentTransformationContainer.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>

PHCosmicsTrkFitter::PHCosmicsTrkFitter(const std::string& name)
  : SubsysReco(name)
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

  // configure alignStates
  m_alignStates.loadNodes(topNode);
  m_alignStates.verbosity(Verbosity());
  m_alignStates.fieldMap(m_fieldMap);

  // detect const field
  std::istringstream stringline(m_fieldMap);
  stringline >> fieldstrength;
  if (!stringline.fail())  // it is a float
  {
    m_ConstField = true;
  }
  auto level = Acts::Logging::FATAL;
  if (Verbosity() > 5)
  {
    level = Acts::Logging::VERBOSE;
  }

  m_fitCfg.fit = ActsTrackFittingAlgorithm::makeKalmanFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField,
      true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", level));

  m_fitCfg.dFit = ActsTrackFittingAlgorithm::makeDirectedKalmanFitterFunction(
      m_tGeometry->geometry().tGeometry,
      m_tGeometry->geometry().magField, true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("DirectedKalman", level));

  MaterialSurfaceSelector selector;
  if (m_directNavigation)
  {
    m_tGeometry->geometry().tGeometry->visitSurfaces(selector, false);
    m_materialSurfaces = selector.surfaces;
  }

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

  _tpccellgeo = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");

  if (m_actsEvaluator)
  {
    m_evaluator = std::make_unique<ActsEvaluator>(m_evalname);
    m_evaluator->Init(topNode);
    m_evaluator->verbosity(Verbosity());
  }

  if (m_seedClusAnalysis)
  {
    m_outfile = new TFile(m_evalname.c_str(), "RECREATE");
    m_tree = new TTree("seedclustree", "Tree with cosmic seeds and their clusters");
    makeBranches();
  }

  if (Verbosity() > 1)
  {
    std::cout << "Finish PHCosmicsTrkFitter Setup" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsTrkFitter::process_event(PHCompositeNode* topNode)
{
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
    {
      logLevel = Acts::Logging::VERBOSE;
    }
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

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsTrkFitter::ResetEvent(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "Reset PHCosmicsTrkFitter" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCosmicsTrkFitter::End(PHCompositeNode* /*topNode*/)
{
  if (m_actsEvaluator)
  {
    m_evaluator->End();
  }
  if (m_seedClusAnalysis)
  {
    m_outfile->cd();
    m_tree->Write();
    m_outfile->Close();
  }
  if (Verbosity() > 0)
  {
    std::cout << "The Acts track fitter succeeded " << m_nGoodFits
              << " times and had " << m_nBadFits
              << " fits return an error" << std::endl;
    std::cout << "There were " << m_nLongSeeds << " long seeds" << std::endl;
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

  for (auto* track : *m_seedMap)
  {
    if (!track)
    {
      continue;
    }

    unsigned int tpcid = track->get_tpc_seed_index();
    unsigned int siid = track->get_silicon_seed_index();

    // get the crossing number
    auto* siseed = m_siliconSeeds->get(siid);
    short crossing = 0;

    auto* tpcseed = m_tpcSeeds->get(tpcid);
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

    if (Verbosity() > 1)
    {
      const auto position = TrackSeedHelper::get_xyz(tpcseed);
      std::cout << " tpc seed position is (x,y,z) = " << position.x() << "  " << position.y() << "  " << position.z() << std::endl;
    }

    ActsTrackFittingAlgorithm::MeasurementContainer measurements;

    SourceLinkVec sourceLinks;

    MakeSourceLinks makeSourceLinks;
    makeSourceLinks.initialize(_tpccellgeo);
    makeSourceLinks.setVerbosity(Verbosity());
    makeSourceLinks.set_pp_mode(false);

    makeSourceLinks.resetTransientTransformMap(
        m_alignmentTransformationMapTransient,
        m_transient_id_set,
        m_tGeometry);

    if (siseed)
    {
      // silicon source links
      sourceLinks = makeSourceLinks.getSourceLinks(
          siseed,
          measurements,
          m_clusterContainer,
          m_tGeometry,
          m_globalPositionWrapper,
          m_alignmentTransformationMapTransient,
          m_transient_id_set,
          crossing);
    }

    // tpc source links
    const auto tpcSourceLinks = makeSourceLinks.getSourceLinks(
        tpcseed,
        measurements,
        m_clusterContainer,
        m_tGeometry,
        m_globalPositionWrapper,
        m_alignmentTransformationMapTransient,
        m_transient_id_set,
        crossing);

    sourceLinks.insert(sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end());

    if (sourceLinks.size() < 5)
    {
      continue;
    }

    Acts::GeometryContext geoContext{m_alignmentTransformationMapTransient};
    // copy transient map for this track into transient geoContext
    m_transient_geocontext = geoContext;

    std::vector<Acts::Vector3> pos;
    std::vector<Acts::Vector3> sorted_positions;
    // get positions from cluster keys
    // TODO: should implement distortions
    TrackSeedHelper::position_map_t positions;
    for (auto key_iter = tpcseed->begin_cluster_keys(); key_iter != tpcseed->end_cluster_keys(); ++key_iter)
    {
      const auto& key(*key_iter);
      positions.emplace(key, m_tGeometry->getGlobalPosition(key, m_clusterContainer->findCluster(key)));
      pos.push_back(positions[key]);
    }
    sorted_positions = pos;

    std::sort(sorted_positions.begin(), sorted_positions.end(), [](const Acts::Vector3& a, const Acts::Vector3& b)
              {
                float aradius = std::sqrt(a.x()*a.x()+a.y()*a.y());
                if(a.y() < 0)
                {
                  aradius *= -1;
                }
                float bradius = std::sqrt(b.x()*b.x()+b.y()*b.y());
                if(b.y() < 0)
                {
                  bradius *= -1;
                }
                return aradius > bradius; });

    TrackSeedHelper::circleFitByTaubin(tpcseed, positions, 0, 58);

    Acts::Vector3 pca = calculatePCA(tpcseed, sorted_positions);

    Acts::Vector3 momentum = calculateMomentum(tpcseed, sorted_positions);

    Acts::Vector3 position = pca * Acts::UnitConstants::cm;

    if (!is_valid(momentum) || !is_valid(position))
    {
      continue;
    }

    int charge = getCharge(tpcseed, sorted_positions);
    
    auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
        position);
    auto actsFourPos = Acts::Vector4(position(0), position(1),
                                     position(2),
                                     10 * Acts::UnitConstants::ns);

    if (sourceLinks.size() > 20)
    {
      m_nLongSeeds++;
    }
    Acts::BoundSquareMatrix cov = setDefaultCovariance();
    if (m_seedClusAnalysis)
    {
      clearVectors();
      m_seed = tpcid;
      m_R = std::abs(1. / tpcseed->get_qOverR());
      m_X0 = tpcseed->get_X0();
      m_Y0 = tpcseed->get_Y0();
      m_Z0 = tpcseed->get_Z0();
      m_slope = tpcseed->get_slope();
      m_pcax = position(0) / Acts::UnitConstants::cm;
      m_pcay = position(1) / Acts::UnitConstants::cm;
      m_pcaz = position(2) / Acts::UnitConstants::cm;
      m_px = momentum(0);
      m_py = momentum(1);
      m_pz = momentum(2);
      
      m_charge = charge;
      fillVectors(tpcseed, siseed);
      m_x.push_back(position.x() / Acts::UnitConstants::cm);
      m_y.push_back(position.y() / Acts::UnitConstants::cm);
      m_z.push_back(position.z() / Acts::UnitConstants::cm);
      m_r.push_back(radius(position.x(), position.y()) / Acts::UnitConstants::cm);
      m_tree->Fill();
    }
    if (m_dumpSeeds)
    {
      SvtxTrack_v4 newTrack;
      newTrack.set_tpc_seed(tpcseed);
      newTrack.set_crossing(crossing);
      newTrack.set_silicon_seed(siseed);

      unsigned int trid = m_trackMap->size();
      newTrack.set_id(trid);
      newTrack.set_px(momentum.x());
      newTrack.set_py(momentum.y());
      newTrack.set_pz(momentum.z());
      newTrack.set_x(position.x());
      newTrack.set_y(position.y());
      newTrack.set_z(position.z());
      newTrack.set_charge(charge);
      m_trackMap->insertWithKey(&newTrack, trid);
      continue;
    }
    //! Reset the track seed with the dummy covariance
    auto seed = ActsTrackFittingAlgorithm::TrackParameters::create(
        m_transient_geocontext,
        pSurface,
        actsFourPos,
        momentum,
        charge / momentum.norm(),
        cov,
        Acts::ParticleHypothesis::muon(),
        100 * Acts::UnitConstants::cm);
    if (!seed.ok())
    {
      std::cout << "Could not create track params, skipping track" << std::endl;
      continue;
    }

    if (Verbosity() > 0)
    {
      std::cout << "Source link size " << sourceLinks.size() << std::endl;
      printTrackSeed(seed.value());
    }

    //! Set host of propagator options for Acts to do e.g. material integration
    auto calibptr = std::make_unique<Calibrator>();
    CalibratorAdapter calibrator{*calibptr, measurements};

    auto magcontext = m_tGeometry->geometry().magFieldContext;
    auto calibcontext = m_tGeometry->geometry().calibContext;
    auto ppPlainOptions = Acts::PropagatorPlainOptions(m_transient_geocontext, magcontext);

    ActsTrackFittingAlgorithm::GeneralFitterOptions
        kfOptions{
            m_transient_geocontext,
            magcontext,
            calibcontext,
            &(*pSurface),
            ppPlainOptions};

    auto trackContainer =
        std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer =
        std::make_shared<Acts::VectorMultiTrajectory>();
    ActsTrackFittingAlgorithm::TrackContainer
        tracks(trackContainer, trackStateContainer);
    auto result = fitTrack(sourceLinks, seed.value(), kfOptions,
                           calibrator, tracks);

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
      m_nGoodFits++;
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

bool PHCosmicsTrkFitter::getTrackFitResult(FitResult& fitOutput,
                                           TrackSeed* seed, SvtxTrack* track,
                                           ActsTrackFittingAlgorithm::TrackContainer& tracks,
                                           const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  /// Make a trajectory state for storage, which conforms to Acts track fit
  /// analysis tool
  auto& outtrack = fitOutput.value();
  std::vector<Acts::TrackIndexType> trackTips;
  trackTips.reserve(1);
  trackTips.emplace_back(outtrack.tipIndex());
  Trajectory::IndexedParameters indexedParams;

  indexedParams.emplace(std::pair{outtrack.tipIndex(),
                                  ActsExamples::TrackParameters{outtrack.referenceSurface().getSharedPtr(),
                                                                outtrack.parameters(), outtrack.covariance(), outtrack.particleHypothesis()}});

  if (Verbosity() > 2)
  {
    std::cout << "Fitted parameters for track" << std::endl;
    std::cout << " position : " << outtrack.referenceSurface().localToGlobal(m_transient_geocontext, Acts::Vector2(outtrack.loc0(), outtrack.loc1()), Acts::Vector3(1, 1, 1)).transpose()

              << std::endl;

    int otcharge = outtrack.qOverP() > 0 ? 1 : -1;
    std::cout << "charge: " << otcharge << std::endl;
    std::cout << " momentum : " << outtrack.momentum().transpose()
              << std::endl;
    std::cout << "For trackTip == " << outtrack.tipIndex() << std::endl;
  }

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
    m_evaluator->evaluateTrackFit(tracks, trackTips, indexedParams, track,
                                  seed, measurements);
  }

  return true;
}

inline ActsTrackFittingAlgorithm::TrackFitterResult PHCosmicsTrkFitter::fitTrack(
    const std::vector<Acts::SourceLink>& sourceLinks,
    const ActsTrackFittingAlgorithm::TrackParameters& seed,
    const ActsTrackFittingAlgorithm::GeneralFitterOptions& kfOptions,
    const CalibratorAdapter& calibrator,
    ActsTrackFittingAlgorithm::TrackContainer& tracks)
{
  return (*m_fitCfg.fit)(sourceLinks, seed, kfOptions, calibrator, tracks);
}

void PHCosmicsTrkFitter::updateSvtxTrack(
    std::vector<Acts::TrackIndexType>& tips,
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
  track->set_x(params.position(m_transient_geocontext)(0) / Acts::UnitConstants::cm);
  track->set_y(params.position(m_transient_geocontext)(1) / Acts::UnitConstants::cm);
  track->set_z(params.position(m_transient_geocontext)(2) / Acts::UnitConstants::cm);

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
                                m_transient_geocontext);
  }

  if (Verbosity() > 2)
  {
    std::cout << " Identify fitted track after updating track states:"
              << std::endl;
    track->identify();
  }

  return;
}

Acts::BoundSquareMatrix PHCosmicsTrkFitter::setDefaultCovariance() const
{
  Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();

  /// Acts cares about the track covariance as it helps the KF
  /// know whether or not to trust the initial track seed or not.
  /// We reset it here to some loose values as it helps Acts improve
  /// the fitting.
  /// If the covariance is too loose, it won't be able to propagate,
  /// but if it is too tight, it will just "believe" the track seed over
  /// the hit data

  // cppcheck-suppress duplicateAssignExpression
  double sigmaD0 = 300 * Acts::UnitConstants::um;
  double sigmaZ0 = 300 * Acts::UnitConstants::um;
  // cppcheck-suppress duplicateAssignExpression
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
      << "position: " << seed.position(m_transient_geocontext).transpose()
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
    auto* node = new PHDataNode<SvtxAlignmentStateMap>(m_alignmentStateMap, "SvtxAlignmentStateMap", "PHObject");
    svtxNode->addNode(node);
  }

  m_alignmentTransformationMapTransient = findNode::getClass<alignmentTransformationContainer>(topNode, "alignmentTransformationContainerTransient");
  if (!m_alignmentTransformationMapTransient)
  {
    std::cout << PHWHERE << "alignmentTransformationContainerTransient not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
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
  // tpc seeds
  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!m_tpcSeeds)
  {
    std::cout << PHWHERE << "TpcTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // silicon seeds
  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!m_siliconSeeds)
  {
    std::cout << PHWHERE << "SiliconTrackSeedContainer not on node tree. Bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // clusters
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE
              << "No trkr cluster container, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsGeometry not on node tree. Exiting."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tracks
  m_seedMap = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if (!m_seedMap)
  {
    std::cout << "No Svtx seed map on node tree. Exiting."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHCosmicsTrkFitter::makeBranches()
{
  m_tree->Branch("seed", &m_seed, "m_seed/I");
  m_tree->Branch("event", &m_event, "m_event/I");
  m_tree->Branch("R", &m_R, "m_R/F");
  m_tree->Branch("X0", &m_X0, "m_X0/F");
  m_tree->Branch("Y0", &m_Y0, "m_Y0/F");
  m_tree->Branch("Z0", &m_Z0, "m_Z0/F");
  m_tree->Branch("slope", &m_slope, "m_slope/F");
  m_tree->Branch("pcax", &m_pcax, "m_pcax/F");
  m_tree->Branch("pcay", &m_pcay, "m_pcay/F");
  m_tree->Branch("pcaz", &m_pcaz, "m_pcaz/F");
  m_tree->Branch("px", &m_px, "m_px/F");
  m_tree->Branch("py", &m_py, "m_py/F");
  m_tree->Branch("pz", &m_pz, "m_pz/F");
  m_tree->Branch("charge", &m_charge, "m_charge/I");
  m_tree->Branch("nmaps", &m_nmaps, "m_nmaps/I");
  m_tree->Branch("nintt", &m_nintt, "m_nintt/I");
  m_tree->Branch("ntpc", &m_ntpc, "m_ntpc/I");
  m_tree->Branch("nmm", &m_nmm, "m_nmm/I");
  m_tree->Branch("locx", &m_locx);
  m_tree->Branch("locy", &m_locy);
  m_tree->Branch("x", &m_x);
  m_tree->Branch("y", &m_y);
  m_tree->Branch("z", &m_z);
  m_tree->Branch("r", &m_r);
  m_tree->Branch("layer", &m_layer);
  m_tree->Branch("phi", &m_phi);
  m_tree->Branch("eta", &m_eta);
  m_tree->Branch("phisize", &m_phisize);
  m_tree->Branch("zsize", &m_zsize);
  m_tree->Branch("ephi", &m_ephi);
  m_tree->Branch("ez", &m_ez);
}
void PHCosmicsTrkFitter::fillVectors(TrackSeed* tpcseed, TrackSeed* siseed)
{
  for (auto* seed : {tpcseed, siseed})
  {
    if (!seed)
    {
      continue;
    }

    for (auto it = seed->begin_cluster_keys(); it != seed->end_cluster_keys();
         ++it)
    {
      auto key = *it;
      auto* cluster = m_clusterContainer->findCluster(key);
      m_locx.push_back(cluster->getLocalX());
      m_locy.push_back(cluster->getLocalY());
      auto glob = m_tGeometry->getGlobalPosition(key, cluster);
      m_x.push_back(glob.x());
      m_y.push_back(glob.y());
      m_z.push_back(glob.z());
      float r = std::sqrt(glob.x() * glob.x() + glob.y() * glob.y());
      if (glob.y() < 0)
      {
        r *= -1;
      }
      m_r.push_back(r);
      TVector3 globt(glob.x(), glob.y(), glob.z());
      m_phi.push_back(globt.Phi());
      m_eta.push_back(globt.Eta());
      m_phisize.push_back(cluster->getPhiSize());
      m_zsize.push_back(cluster->getZSize());
      auto para_errors =
          ClusterErrorPara::get_clusterv5_modified_error(cluster, r, key);

      m_ephi.push_back(std::sqrt(para_errors.first));
      m_ez.push_back(std::sqrt(para_errors.second));
    }
  }
}
void PHCosmicsTrkFitter::clearVectors()
{
  m_locx.clear();
  m_locy.clear();
  m_x.clear();
  m_y.clear();
  m_z.clear();
  m_layer.clear();
  m_r.clear();
  m_phi.clear();
  m_eta.clear();
  m_phisize.clear();
  m_zsize.clear();
  m_ephi.clear();
  m_ez.clear();
}

int PHCosmicsTrkFitter::getCharge(TrackSeed* tpcseed,
                                  const std::vector<Acts::Vector3>& sorted_positions)
{
  Acts::GeometryContext transient_geocontext{m_alignmentTransformationMapTransient};

  std::vector<float> tpcparams{(float) std::abs(1. / tpcseed->get_qOverR()),
                               tpcseed->get_X0(),
                               tpcseed->get_Y0(),
                               tpcseed->get_slope(),
                               tpcseed->get_Z0()};

  float phi0 = std::atan2(sorted_positions[0].y() - tpcparams[2], sorted_positions[0].x() - tpcparams[1]);
  int posphi = 0;
  int negphi = 0;
  // just take the first 4 outermost clusters as a test to determine the bend angle
  // from the outermost radial cluster
  for (size_t i = 1; i < 5; i++)
  {
    auto cluspos = sorted_positions[i];

    float phi = std::atan2(cluspos.y() - tpcparams[2], cluspos.x() - tpcparams[1]);
    if (phi > phi0)
    {
      posphi++;
    }
    else
    {
      negphi++;
    }
  }
  int charge = posphi > negphi ? -1 : 1;
  if (Verbosity() > 2)
  {
    std::cout << "charge is " << charge << std::endl;
  }

  return charge;
}

Acts::Vector3 PHCosmicsTrkFitter::calculatePCA(TrackSeed* seed, const std::vector<Acts::Vector3>& sorted_positions) const
{
  float tpcR = fabs(1. / seed->get_qOverR());
  float tpcx = seed->get_X0();
  float tpcy = seed->get_Y0();

  // calculate the pcaxy for the seed wrt a line surface located at (0,m_vertexRadius) in x-y plane
  float dx = -tpcx;
  float dy = m_vertexRadius - tpcy;
  float dist = std::sqrt(dx * dx + dy * dy);
  float pcax = tpcx + tpcR * (dx / dist);
  float pcay = tpcy + tpcR * (dy / dist);

  auto arcLength = [&](float x, float y)
  {
    float angle = std::atan2(y - tpcy, x - tpcx);
    return tpcR * angle;
  };

  float sum_s = 0;
  float sum_z = 0;
  float sum_ss = 0;
  float sum_sz = 0;
  int n = sorted_positions.size();
  // Compute the arc-length parameter for each cluster, then fit to a line
  // Fit z = a + b*s using simple linear regression
  for (const auto& p : sorted_positions)
  {
    float s = arcLength(p.x(), p.y());
    sum_s += s;
    sum_z += p.z();
    sum_ss += s * s;
    sum_sz += s * p.z();
  }

  float denom = n * sum_ss - sum_s * sum_s;
  float b = (n * sum_sz - sum_s * sum_z) / denom;
  float a = (sum_z - b * sum_s) / n;

  // Then evaluate at the arc length of the PCA to get the z position of the PCA
  float s_ca = arcLength(pcax, pcay);
  float z_ca = a + b * s_ca;

  return Acts::Vector3(pcax, pcay, z_ca);
}

Acts::Vector3 PHCosmicsTrkFitter::calculateMomentum(TrackSeed* tpcseed, const std::vector<Acts::Vector3>& sorted_positions)
{
  // now calculate the momentum vector
  const auto intersect = TrackFitUtils::circle_circle_intersection(m_vertexRadius, std::abs(1. / tpcseed->get_qOverR()), tpcseed->get_X0(), tpcseed->get_Y0());
  float intx;
  float inty;

  if (std::get<1>(intersect) > std::get<3>(intersect))
  {
    intx = std::get<0>(intersect);
    inty = std::get<1>(intersect);
  }
  else
  {
    intx = std::get<2>(intersect);
    inty = std::get<3>(intersect);
  }
  if (Verbosity() > 2)
  {
    std::cout << "XY intersection options " << std::get<0>(intersect) << ", " << std::get<1>(intersect) << " and " << std::get<2>(intersect) << ", " << std::get<3>(intersect) << std::endl;
  }

  TrackFitUtils::position_vector_t xypoints;
  TrackFitUtils::position_vector_t rzpoints;
  for (const auto& p : sorted_positions)
  {
    float clusr = radius(p.x(), p.y());
    if (p.y() < 0)
    {
      clusr *= -1;
    }

    // exclude silicon and tpot clusters for now
    if (std::abs(clusr) > 80 || std::abs(clusr) < 30)
    {
      continue;
    }
    xypoints.emplace_back(p.x(), p.y());
    rzpoints.emplace_back(p.z(), clusr);
  }

  auto rzparams = TrackFitUtils::line_fit(rzpoints);
  float fulllineslope = std::get<0>(rzparams);

  float slope = tpcseed->get_slope();
  float intz = m_vertexRadius * slope + tpcseed->get_Z0();

  Acts::Vector3 inter(intx, inty, intz);

  std::vector<float> tpcparams{(float) std::abs(1. / tpcseed->get_qOverR()),
                               tpcseed->get_X0(),
                               tpcseed->get_Y0(),
                               tpcseed->get_slope(),
                               tpcseed->get_Z0()};
  auto tangent = TrackFitUtils::get_helix_tangent(tpcparams,
                                                  inter);

  auto tan = tangent.second;

  float p;
  if (m_ConstField)
  {
    p = std::cosh(tpcseed->get_eta()) * fabs(1. / tpcseed->get_qOverR()) * (0.3 / 100) * fieldstrength;
  }
  else
  {
    p = tpcseed->get_p();
  }

  tan *= p;

  //! if we got the opposite seed then z will be backwards, so we take the
  //! value of tan.z() multiplied by the sign of the slope determined for
  //! the full cosmic track
  //! same with px/py since a single cosmic produces two seeds that bend
  //! in opposite directions
  float theta = std::atan(fulllineslope);
  /// Normalize to 0<theta<pi
  if (theta < 0)
  {
    theta += M_PI;
  }
  float pz = std::cos(theta) * p;
  if (fulllineslope < 0)
  {
    pz = std::abs(pz);
  }
  else
  {
    pz = std::abs(pz) * -1;
  }
  Acts::Vector3 momentum = Acts::Vector3::Zero();

  if (!m_zeroField)
  {
    momentum.x() = tan.x() * -1;
    momentum.y() = tan.y() * -1;
  }
  else
  {
    auto xyparams = TrackFitUtils::line_fit(xypoints);
    float fulllineslopexy = std::get<0>(xyparams);
    if (fulllineslopexy < 0)
    {
      momentum.x() = fabs(tan.x());
    }
    else
    {
      momentum.x() = fabs(tan.x()) * -1;
    }
    momentum.y() = fabs(tan.y()) * -1;
  }

  momentum.z() = pz;

  return momentum;
}