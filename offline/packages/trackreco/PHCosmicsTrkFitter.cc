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
#include <TVector3.h>

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

  _tpccellgeo = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

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

  m_event++;
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

    if (Verbosity() > 1)
    {
      std::cout << " tpc seed position is (x,y,z) = " << tpcseed->get_x() << "  " << tpcseed->get_y() << "  " << tpcseed->get_z() << std::endl;
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
      sourceLinks = makeSourceLinks.getSourceLinks(
          siseed,
          measurements,
          m_clusterContainer,
          m_tGeometry,
          _dcc_static, _dcc_average, _dcc_fluctuation,
          m_alignmentTransformationMapTransient,
          m_transient_id_set,
          crossing);
    }
    const auto tpcSourceLinks = makeSourceLinks.getSourceLinks(
        tpcseed,
        measurements,
        m_clusterContainer,
        m_tGeometry,
        _dcc_static, _dcc_average, _dcc_fluctuation,
        m_alignmentTransformationMapTransient,
        m_transient_id_set,
        crossing);

    sourceLinks.insert(sourceLinks.end(), tpcSourceLinks.begin(), tpcSourceLinks.end());

    if (sourceLinks.size() < 5)
    {
      continue;
    }
    int charge = 0;
    float cosmicslope = 0;

    getCharge(tpcseed, charge, cosmicslope);

    // copy transient map for this track into transient geoContext
    m_transient_geocontext = m_alignmentTransformationMapTransient;

    tpcseed->circleFitByTaubin(m_clusterContainer, m_tGeometry, 0, 58);

    float tpcR = fabs(1. / tpcseed->get_qOverR());
    float tpcx = tpcseed->get_X0();
    float tpcy = tpcseed->get_Y0();

    const auto intersect =
        TrackFitUtils::circle_circle_intersection(m_vertexRadius,
                                                  tpcR, tpcx, tpcy);
    float intx, inty;

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
    std::vector<TrkrDefs::cluskey> keys;
    std::vector<Acts::Vector3> clusPos;
    std::copy(tpcseed->begin_cluster_keys(), tpcseed->end_cluster_keys(), std::back_inserter(keys));
    TrackFitUtils::getTrackletClusters(m_tGeometry, m_clusterContainer,
                                       clusPos, keys);
    TrackFitUtils::position_vector_t xypoints, rzpoints;
    for (auto& pos : clusPos)
    {
      float clusr = radius(pos.x(), pos.y());
      if (pos.y() < 0)
      {
        clusr *= -1;
      }

      // exclude silicon and tpot clusters for now
      if (fabs(clusr) > 80 || fabs(clusr) < 30)
      {
        continue;
      }
      xypoints.push_back(std::make_pair(pos.x(), pos.y()));
      rzpoints.push_back(std::make_pair(pos.z(), clusr));
    }

    auto rzparams = TrackFitUtils::line_fit(rzpoints);
    float fulllineintz = std::get<1>(rzparams);
    float fulllineslope = std::get<0>(rzparams);

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
    float theta = std::atan(fulllineslope);
    /// Normalize to 0<theta<pi
    if (theta < 0)
    {
      theta += M_PI;
    }
    float pz = std::cos(theta) * p;
    if (fulllineslope < 0)
    {
      pz = fabs(pz);
    }
    else
    {
      pz = fabs(pz) * -1;
    }
    Acts::Vector3 momentum = Acts::Vector3::Zero();

    if (!m_zeroField)
    {
      momentum.x() = charge < 0 ? tan.x() : tan.x() * -1;
      momentum.y() = charge < 0 ? tan.y() : tan.y() * -1;
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
    Acts::Vector3 position(pca.x(), pca.y(),
                           (m_vertexRadius - fulllineintz) / fulllineslope);

    position *= Acts::UnitConstants::cm;
    if (!is_valid(momentum))
    {
      continue;
    }

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
      m_R = tpcR;
      m_X0 = tpcx;
      m_Y0 = tpcy;
      m_Z0 = fulllineintz;
      m_slope = fulllineslope;
      m_pcax = position(0);
      m_pcay = position(1);
      m_pcaz = position(2);
      m_px = momentum(0);
      m_py = momentum(1);
      m_pz = momentum(2);
      m_charge = charge;
      fillVectors(siseed, tpcseed);
      m_tree->Fill();
    }
    //! Reset the track seed with the dummy covariance
    auto seed = ActsTrackFittingAlgorithm::TrackParameters::create(
        pSurface,
        m_transient_geocontext,
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
    Acts::PropagatorPlainOptions ppPlainOptions;

    auto calibptr = std::make_unique<Calibrator>();
    CalibratorAdapter calibrator{*calibptr, measurements};

    auto magcontext = m_tGeometry->geometry().magFieldContext;
    auto calibcontext = m_tGeometry->geometry().calibContext;

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
  std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
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
  for (auto seed : {tpcseed, siseed})
  {
    if (!seed)
    {
      continue;
    }

    for (auto it = seed->begin_cluster_keys(); it != seed->end_cluster_keys();
         ++it)
    {
      auto key = *it;
      auto cluster = m_clusterContainer->findCluster(key);
      m_locx.push_back(cluster->getLocalX());
      float ly = cluster->getLocalY();
      if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::tpcId)
      {
        double drift_velocity = m_tGeometry->get_drift_velocity();
        double zdriftlength = cluster->getLocalY() * drift_velocity;
        double surfCenterZ = 52.89;                // 52.89 is where G4 thinks the surface center is
        double zloc = surfCenterZ - zdriftlength;  // converts z drift length to local z position in the TPC in north
        unsigned int side = TpcDefs::getSide(key);
        if (side == 0) zloc = -zloc;
        ly = zloc * 10;
      }
      m_locy.push_back(ly);
      auto glob = m_tGeometry->getGlobalPosition(key, cluster);
      m_x.push_back(glob.x());
      m_y.push_back(glob.y());
      m_z.push_back(glob.z());
      float r = std::sqrt(glob.x() * glob.x() + glob.y() * glob.y());
      m_r.push_back(r);
      TVector3 globt(glob.x(), glob.y(), glob.z());
      m_phi.push_back(globt.Phi());
      m_eta.push_back(globt.Eta());
      m_phisize.push_back(cluster->getPhiSize());
      m_zsize.push_back(cluster->getZSize());
      auto para_errors =
          m_clusErrPara.get_clusterv5_modified_error(cluster, r, key);

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

void PHCosmicsTrkFitter::getCharge(
    TrackSeed* track,
    // TrkrClusterContainer*  clusterContainer,
    // ActsGeometry* tGeometry,
    // alignmentTransformationContainer* transformMapTransient,
    // float vertexRadius,
    int& charge,
    float& cosmicslope)
{
  Acts::GeometryContext transient_geocontext;
  transient_geocontext = m_alignmentTransformationMapTransient;  // set local/global transforms to distortion corrected ones for this track

  std::vector<Acts::Vector3> global_vec;

  for (auto clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
  {
    auto key = *clusIter;
    auto cluster = m_clusterContainer->findCluster(key);
    if (!cluster)
    {
      std::cout << "MakeSourceLinks::getCharge: Failed to get cluster with key " << key << " for track seed" << std::endl;
      continue;
    }

    auto surf = m_tGeometry->maps().getSurface(key, cluster);
    if (!surf)
    {
      continue;
    }

    // get cluster global positions
    Acts::Vector2 local = m_tGeometry->getLocalCoords(key, cluster);  // converts TPC time to z
    Acts::Vector3 glob = surf->localToGlobal(transient_geocontext,
                                             local * Acts::UnitConstants::cm,
                                             Acts::Vector3(1, 1, 1));
    glob /= Acts::UnitConstants::cm;

    global_vec.push_back(glob);
  }

  Acts::Vector3 globalMostOuter;
  Acts::Vector3 globalSecondMostOuter(0, 999999, 0);
  float largestR = 0;
  // loop over global positions
  for (int i = 0; i < global_vec.size(); ++i)
  {
    Acts::Vector3 global = global_vec[i];
    // float r = std::sqrt(square(global.x()) + square(global.y()));
    float r = radius(global.x(), global.y());

    /// use the top hemisphere to determine the charge
    if (r > largestR && global.y() > 0)
    {
      globalMostOuter = global_vec[i];
      largestR = r;
    }
  }

  //! find the closest cluster to the outermost cluster
  float maxdr = std::numeric_limits<float>::max();
  for (int i = 0; i < global_vec.size(); i++)
  {
    if (global_vec[i].y() < 0) continue;

    float dr = std::sqrt(square(globalMostOuter.x()) + square(globalMostOuter.y())) - std::sqrt(square(global_vec[i].x()) + square(global_vec[i].y()));
    //! Place a dr cut to get maximum bend due to TPC clusters having
    //! larger fluctuations
    if (dr < maxdr && dr > 10)
    {
      maxdr = dr;
      globalSecondMostOuter = global_vec[i];
    }
  }

  //! we have to calculate phi WRT the vertex position outside the detector,
  //! not at (0,0)
  Acts::Vector3 vertex(0, m_vertexRadius, 0);
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

  return;
}
