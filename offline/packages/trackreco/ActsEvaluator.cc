#include "ActsEvaluator.h"

/// General fun4all and subsysreco includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <g4main/PHG4VtxPoint.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/MultiTrajectoryHelpers.hpp>

#include <TFile.h>
#include <TTree.h>

ActsEvaluator::ActsEvaluator(const std::string& name)
  : m_filename(name)
{
}

ActsEvaluator::~ActsEvaluator()
{
}

void ActsEvaluator::Init(PHCompositeNode* topNode)
{
  if (m_verbosity > 1)
  {
    std::cout << "Starting ActsEvaluator::Init" << std::endl;
  }

  if (getNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
  {
    std::cout << "Error: Nodes not available in ActsEvaluator::Init"
              << std::endl;
    return;
  }

  initializeTree();

  if (m_verbosity > 1)
  {
    std::cout << "Finished ActsEvaluator::Init" << std::endl;
  }
}
void ActsEvaluator::next_event(PHCompositeNode* topNode)
{
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack = new SvtxEvalStack(topNode);
  }

  m_svtxEvalStack->next_event(topNode);
  
  m_eventNr++;
}
void ActsEvaluator::process_track(const ActsTrackFittingAlgorithm::TrackContainer& tracks,
				  std::vector<Acts::MultiTrajectoryTraits::IndexType>& trackTips,
				  Trajectory::IndexedParameters& paramsMap,
                                  SvtxTrack* track,
                                  const TrackSeed* seed,
                                  const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  if (m_verbosity > 1)
  {
    std::cout << "Starting ActsEvaluator at event " << m_eventNr
              << std::endl;
  }

  evaluateTrackFit(tracks,trackTips,paramsMap, track, seed, measurements);


  if (m_verbosity > 1)
  {
    std::cout << "Finished ActsEvaluator::process_event" << std::endl;
  }
}

void ActsEvaluator::evaluateTrackFit(const ActsTrackFittingAlgorithm::TrackContainer& tracks,
				     std::vector<Acts::MultiTrajectoryTraits::IndexType>& trackTips,
				     Trajectory::IndexedParameters& paramsMap,
                                     SvtxTrack* track,
                                     const TrackSeed* seed,
                                     const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  if (m_verbosity > 5)
    {
      std::cout << "Evaluating Acts track fits" << std::endl;
    }
  /// Skip failed fits
  if (trackTips.empty())
  {
    if (m_verbosity > 1)
    {
      std::cout << "TrackTips empty in ActsEvaluator" << std::endl;
    }
    return;
  }

  SvtxTrackEval* trackeval = m_svtxEvalStack->get_track_eval();

  int iTrack = track->get_id();
  int iTraj = iTrack;
  if (m_verbosity > 2)
  {
    std::cout << "Starting trajectory with trackKey " << track->get_id()
              << " and corresponding to tpc track seed "
              << seed->get_tpc_seed_index() << std::endl;
  }

  const auto& mj = tracks.trackStateContainer();
  const auto& trackTip = trackTips.front();
  m_trajNr = iTraj;

  if (m_verbosity > 2)
  {
    std::cout << "beginning trackTip " << trackTip
              << " corresponding to track " << iTrack
              << " in trajectroy " << iTraj << std::endl;
  }

  if (m_verbosity > 2)
  {
    std::cout << "Evaluating track key " << iTrack
              << " for track tip " << trackTip << std::endl;
  }

  PHG4Particle* g4particle = trackeval->max_truth_particle_by_nclusters(track);

  if (m_verbosity > 1)
  {
    std::cout << "Analyzing SvtxTrack " << iTrack << std::endl;

    std::cout << "TruthParticle : " << g4particle->get_px()
              << ", " << g4particle->get_py() << ", "
              << g4particle->get_pz() << ", " << g4particle->get_e()
              << std::endl;
  }

  m_trackNr = iTrack;

  auto trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
  
  const auto& params = paramsMap.find(trackTip)->second;

  if (m_verbosity > 1)
  {
   
      std::cout << "Fitted params : "
                << params.position(m_tGeometry->geometry().getGeoContext())
                << std::endl
                << params.momentum()
                << std::endl;
      std::cout << "Track has " << trajState.nMeasurements
                << " measurements and " << trajState.nHoles
                << " holes and " << trajState.nOutliers
                << " outliers and " << trajState.nStates
                << " states " << std::endl;
    
  }

  m_nMeasurements = trajState.nMeasurements;
  m_nStates = trajState.nStates;
  m_nOutliers = trajState.nOutliers;
  m_nSharedHits = trajState.nSharedHits;
  m_nHoles = trajState.nHoles;
  m_chi2_fit = trajState.chi2Sum;
  m_ndf_fit = trajState.NDF;
  m_quality = track->get_quality();

  fillG4Particle(g4particle);
  fillProtoTrack(seed);
  fillFittedTrackParams(paramsMap, trackTip);
  visitTrackStates(mj, trackTip, measurements);

  m_trackTree->Fill();

  /// Start fresh for the next track
  clearTrackVariables();
  if (m_verbosity > 1)
  {
    std::cout << "Finished track " << iTrack << std::endl;
  }

  if (m_verbosity > 1)
  {
    std::cout << "Analyzed " << iTrack << " tracks in trajectory number "
              << iTraj << std::endl;
  }

  if (m_verbosity > 5)
  {
    std::cout << "Finished evaluating track fits" << std::endl;
  }
  return;
}

void ActsEvaluator::End()
{
  m_trackFile->cd();
  m_trackTree->Write();
  m_trackFile->Close();
}

void ActsEvaluator::visitTrackStates(const Acts::ConstVectorMultiTrajectory& traj,
                                     const size_t& trackTip,
                                     const ActsTrackFittingAlgorithm::MeasurementContainer& measurements)
{
  if (m_verbosity > 2)
  {
    std::cout << "Begin visit track states" << std::endl;
  }

  traj.visitBackwards(trackTip, [&](const auto& state)
                      {
    /// Only fill the track states with non-outlier measurement
    auto typeFlags = state.typeFlags();
    if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag))
    {
      return true;
    }

    const auto& surface = state.referenceSurface();

    /// Get the geometry ID
    auto geoID = surface.geometryId();
    m_volumeID.push_back(geoID.volume());
    m_layerID.push_back(geoID.layer());
    m_moduleID.push_back(geoID.sensitive());

    if(m_verbosity > 3)
      {
      std::cout << "Cluster volume : layer : sensitive " << geoID.volume()
		<< " : " << geoID.layer() << " : " 
		<< geoID.sensitive() << std::endl;
      }
    auto sourceLink = state.getUncalibratedSourceLink().template get<ActsSourceLink>();
    const auto& cluskey = sourceLink.cluskey();
        
    Acts::Vector2 local = Acts::Vector2::Zero();
   
    /// get the local measurement that acts used
    std::visit([&](const auto& meas) {
	local(0) = meas.parameters()[0];
	local(1) = meas.parameters()[1];
      }, measurements[sourceLink.index()]);
    
    /// Get global position
    /// This is an arbitrary vector. Doesn't matter in coordinate transformation
    /// in Acts code
    Acts::Vector3 mom(1, 1, 1);
    Acts::Vector3 global = surface.localToGlobal(m_tGeometry->geometry().getGeoContext(),
						 local, mom);

    /// Get measurement covariance
    /// These are the same
    //auto slcov = sourceLink.covariance();
    auto cov = state.effectiveCalibratedCovariance();

    m_lx_hit.push_back(local.x());
    m_ly_hit.push_back(local.y());
    m_x_hit.push_back(global.x());
    m_y_hit.push_back(global.y());
    m_z_hit.push_back(global.z());

    /// Get the truth hit corresponding to this trackState
    /// We go backwards from hitID -> TrkrDefs::cluskey to g4hit with
    /// the map created in PHActsSourceLinks
    float gt = -9999;
    Acts::Vector3 globalTruthPos = getGlobalTruthHit(cluskey, gt);
    float gx = globalTruthPos(0);
    float gy = globalTruthPos(1);
    float gz = globalTruthPos(2);

    /// Get local truth position
    const float r = sqrt(gx * gx + gy * gy + gz * gz);
    Acts::Vector3 globalTruthUnitDir(gx / r, gy / r, gz / r);

    auto vecResult = surface.globalToLocal(
        m_tGeometry->geometry().getGeoContext(),
        globalTruthPos,
        globalTruthUnitDir);

    /// Push the truth hit info
    m_t_x.push_back(gx);
    m_t_y.push_back(gy);
    m_t_z.push_back(gz);
    m_t_r.push_back(sqrt(gx * gx + gy * gy));
    m_t_dx.push_back(gx / r);
    m_t_dy.push_back(gy / r);
    m_t_dz.push_back(gz / r);

    /// Get the truth track parameter at this track State
    float truthLOC0 = 0;
    float truthLOC1 = 0;
    float truthPHI = 0;
    float truthTHETA = 0;
    float truthQOP = 0;
    float truthTIME = 0;
    float momentum = sqrt(m_t_px * m_t_px +
                          m_t_py * m_t_py +
                          m_t_pz * m_t_pz);
    
    if(vecResult.ok())
      {
	Acts::Vector2 truthLocVec = vecResult.value();
	truthLOC0 = truthLocVec.x();
	truthLOC1 = truthLocVec.y();
      }
    else
      {
	truthLOC0 = -9999.;
	truthLOC1 = -9999.;
      }

    truthPHI = phi(globalTruthUnitDir);
    truthTHETA = theta(globalTruthUnitDir);
    truthQOP =
        m_t_charge / momentum;
    truthTIME = gt;

    m_t_eLOC0.push_back(truthLOC0);
    m_t_eLOC1.push_back(truthLOC1);
    m_t_ePHI.push_back(truthPHI);
    m_t_eTHETA.push_back(truthTHETA);
    m_t_eQOP.push_back(truthQOP);
    m_t_eT.push_back(truthTIME);

    /// Get the predicted parameter for this state
    bool predicted = false;
    if (state.hasPredicted())
    {
      predicted = true;
      m_nPredicted++;
      
      auto parameters = state.predicted();
      auto covariance = state.predictedCovariance();

      /// Local hit residual info
      auto H = state.effectiveProjector();
      auto resCov = cov + H * covariance * H.transpose();
      auto residual = state.effectiveCalibrated() - H * parameters;
      m_res_x_hit.push_back(residual(Acts::eBoundLoc0));
      m_res_y_hit.push_back(residual(Acts::eBoundLoc1));
      m_err_x_hit.push_back(
          sqrt(resCov(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_err_y_hit.push_back(
          sqrt(resCov(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_pull_x_hit.push_back(
          residual(Acts::eBoundLoc0) /
          sqrt(resCov(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_pull_y_hit.push_back(
          residual(Acts::eBoundLoc1) /
          sqrt(resCov(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_dim_hit.push_back(state.calibratedSize());

      /// Predicted parameter
      m_eLOC0_prt.push_back(parameters[Acts::eBoundLoc0]);
      m_eLOC1_prt.push_back(parameters[Acts::eBoundLoc1]);
      m_ePHI_prt.push_back(parameters[Acts::eBoundPhi]);
      m_eTHETA_prt.push_back(parameters[Acts::eBoundTheta]);
      m_eQOP_prt.push_back(parameters[Acts::eBoundQOverP]);
      m_eT_prt.push_back(parameters[Acts::eBoundTime]);

      /// Predicted residual
      m_res_eLOC0_prt.push_back(parameters[Acts::eBoundLoc0] -
                                truthLOC0);
      m_res_eLOC1_prt.push_back(parameters[Acts::eBoundLoc1] -
                                truthLOC1);
      m_res_ePHI_prt.push_back(parameters[Acts::eBoundPhi] -
                               truthPHI);
      m_res_eTHETA_prt.push_back(
          parameters[Acts::eBoundTheta] - truthTHETA);
      m_res_eQOP_prt.push_back(parameters[Acts::eBoundQOverP] -
                               truthQOP);
      m_res_eT_prt.push_back(parameters[Acts::eBoundTime] -
                             truthTIME);

      /// Predicted parameter Uncertainties
      m_err_eLOC0_prt.push_back(
          sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_err_eLOC1_prt.push_back(
          sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_err_ePHI_prt.push_back(
          sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
      m_err_eTHETA_prt.push_back(
          sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
      m_err_eQOP_prt.push_back(
          sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
      m_err_eT_prt.push_back(
          sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

      /// Predicted parameter pulls
      m_pull_eLOC0_prt.push_back(
          (parameters[Acts::eBoundLoc0] - truthLOC0) /
          sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_pull_eLOC1_prt.push_back(
          (parameters[Acts::eBoundLoc1] - truthLOC1) /
          sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_pull_ePHI_prt.push_back(
          (parameters[Acts::eBoundPhi] - truthPHI) /
          sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
      m_pull_eTHETA_prt.push_back(
          (parameters[Acts::eBoundTheta] - truthTHETA) /
          sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
      m_pull_eQOP_prt.push_back(
          (parameters[Acts::eBoundQOverP] - truthQOP) /
          sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
      m_pull_eT_prt.push_back(
          (parameters[Acts::eBoundTime] - truthTIME) /
          sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

      Acts::FreeVector freeParams = 
	Acts::detail::transformBoundToFreeParameters(state.referenceSurface(), 
						     m_tGeometry->geometry().getGeoContext(),
						     parameters);

      m_x_prt.push_back(freeParams[Acts::eFreePos0]);
      m_y_prt.push_back(freeParams[Acts::eFreePos1]);
      m_z_prt.push_back(freeParams[Acts::eFreePos2]);
      auto p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
      m_px_prt.push_back(p * freeParams[Acts::eFreeDir0]);
      m_py_prt.push_back(p * freeParams[Acts::eFreeDir1]);
      m_pz_prt.push_back(p * freeParams[Acts::eFreeDir2]);
      m_pT_prt.push_back(p * std::hypot(freeParams[Acts::eFreeDir0],
					freeParams[Acts::eFreeDir1]));
      m_eta_prt.push_back(
         Acts::VectorHelpers::eta(freeParams.segment<3>(Acts::eFreeDir0)));
    }
    else
    {
      /// Push bad values if no predicted parameter
      m_res_x_hit.push_back(-9999);
      m_res_y_hit.push_back(-9999);
      m_err_x_hit.push_back(-9999);
      m_err_y_hit.push_back(-9999);
      m_pull_x_hit.push_back(-9999);
      m_pull_y_hit.push_back(-9999);
      m_dim_hit.push_back(-9999);
      m_eLOC0_prt.push_back(-9999);
      m_eLOC1_prt.push_back(-9999);
      m_ePHI_prt.push_back(-9999);
      m_eTHETA_prt.push_back(-9999);
      m_eQOP_prt.push_back(-9999);
      m_eT_prt.push_back(-9999);
      m_res_eLOC0_prt.push_back(-9999);
      m_res_eLOC1_prt.push_back(-9999);
      m_res_ePHI_prt.push_back(-9999);
      m_res_eTHETA_prt.push_back(-9999);
      m_res_eQOP_prt.push_back(-9999);
      m_res_eT_prt.push_back(-9999);
      m_err_eLOC0_prt.push_back(-9999);
      m_err_eLOC1_prt.push_back(-9999);
      m_err_ePHI_prt.push_back(-9999);
      m_err_eTHETA_prt.push_back(-9999);
      m_err_eQOP_prt.push_back(-9999);
      m_err_eT_prt.push_back(-9999);
      m_pull_eLOC0_prt.push_back(-9999);
      m_pull_eLOC1_prt.push_back(-9999);
      m_pull_ePHI_prt.push_back(-9999);
      m_pull_eTHETA_prt.push_back(-9999);
      m_pull_eQOP_prt.push_back(-9999);
      m_pull_eT_prt.push_back(-9999);
      m_x_prt.push_back(-9999);
      m_y_prt.push_back(-9999);
      m_z_prt.push_back(-9999);
      m_px_prt.push_back(-9999);
      m_py_prt.push_back(-9999);
      m_pz_prt.push_back(-9999);
      m_pT_prt.push_back(-9999);
      m_eta_prt.push_back(-9999);
    }

    bool filtered = false;
    if (state.hasFiltered())
    {
      filtered = true;
      m_nFiltered++;

      auto parameter = state.filtered();
      auto covariance = state.filteredCovariance();

      m_eLOC0_flt.push_back(parameter[Acts::eBoundLoc0]);
      m_eLOC1_flt.push_back(parameter[Acts::eBoundLoc1]);
      m_ePHI_flt.push_back(parameter[Acts::eBoundPhi]);
      m_eTHETA_flt.push_back(parameter[Acts::eBoundTheta]);
      m_eQOP_flt.push_back(parameter[Acts::eBoundQOverP]);
      m_eT_flt.push_back(parameter[Acts::eBoundTime]);

      m_res_eLOC0_flt.push_back(parameter[Acts::eBoundLoc0] -
                                truthLOC0);
      m_res_eLOC1_flt.push_back(parameter[Acts::eBoundLoc1] -
                                truthLOC1);
      m_res_ePHI_flt.push_back(parameter[Acts::eBoundPhi] -
                               truthPHI);
      m_res_eTHETA_flt.push_back(parameter[Acts::eBoundTheta] -
                                 truthTHETA);
      m_res_eQOP_flt.push_back(parameter[Acts::eBoundQOverP] -
                               truthQOP);
      m_res_eT_flt.push_back(parameter[Acts::eBoundTime] -
                             truthTIME);

      /// Filtered parameter uncertainties
      m_err_eLOC0_flt.push_back(
          sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_err_eLOC1_flt.push_back(
          sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_err_ePHI_flt.push_back(
          sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
      m_err_eTHETA_flt.push_back(
          sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
      m_err_eQOP_flt.push_back(
          sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
      m_err_eT_flt.push_back(
          sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

      /// Filtered parameter pulls
      m_pull_eLOC0_flt.push_back(
          (parameter[Acts::eBoundLoc0] - truthLOC0) /
          sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_pull_eLOC1_flt.push_back(
          (parameter[Acts::eBoundLoc1] - truthLOC1) /
          sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_pull_ePHI_flt.push_back(
          (parameter[Acts::eBoundPhi] - truthPHI) /
          sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
      m_pull_eTHETA_flt.push_back(
          (parameter[Acts::eBoundTheta] - truthTHETA) /
          sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
      m_pull_eQOP_flt.push_back(
          (parameter[Acts::eBoundQOverP] - truthQOP) /
          sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
      m_pull_eT_flt.push_back(
          (parameter[Acts::eBoundTime] - truthTIME) /
          sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

      Acts::FreeVector freeparams = Acts::detail::transformBoundToFreeParameters(surface, m_tGeometry->geometry().getGeoContext(), parameter);

      /// Other filtered parameter info
      m_x_flt.push_back(freeparams[Acts::eFreePos0]);
      m_y_flt.push_back(freeparams[Acts::eFreePos1]);
      m_z_flt.push_back(freeparams[Acts::eFreePos2]);
      auto p = std::abs(1/freeparams[Acts::eFreeQOverP]);
      m_px_flt.push_back(p * freeparams[Acts::eFreeDir0]);
      m_py_flt.push_back(p * freeparams[Acts::eFreeDir1]);
      m_pz_flt.push_back(p * freeparams[Acts::eFreeDir2]);
      m_pT_flt.push_back(sqrt(p * std::hypot(freeparams[Acts::eFreeDir0],
					     freeparams[Acts::eFreeDir1])));
      m_eta_flt.push_back(eta(freeparams.segment<3>(Acts::eFreeDir0)));
      m_chi2.push_back(state.chi2());
      
    }
    else
    {
      /// Push bad values if no filtered parameter
      m_eLOC0_flt.push_back(-9999);
      m_eLOC1_flt.push_back(-9999);
      m_ePHI_flt.push_back(-9999);
      m_eTHETA_flt.push_back(-9999);
      m_eQOP_flt.push_back(-9999);
      m_eT_flt.push_back(-9999);
      m_res_eLOC0_flt.push_back(-9999);
      m_res_eLOC1_flt.push_back(-9999);
      m_res_ePHI_flt.push_back(-9999);
      m_res_eTHETA_flt.push_back(-9999);
      m_res_eQOP_flt.push_back(-9999);
      m_res_eT_flt.push_back(-9999);
      m_err_eLOC0_flt.push_back(-9999);
      m_err_eLOC1_flt.push_back(-9999);
      m_err_ePHI_flt.push_back(-9999);
      m_err_eTHETA_flt.push_back(-9999);
      m_err_eQOP_flt.push_back(-9999);
      m_err_eT_flt.push_back(-9999);
      m_pull_eLOC0_flt.push_back(-9999);
      m_pull_eLOC1_flt.push_back(-9999);
      m_pull_ePHI_flt.push_back(-9999);
      m_pull_eTHETA_flt.push_back(-9999);
      m_pull_eQOP_flt.push_back(-9999);
      m_pull_eT_flt.push_back(-9999);
      m_x_flt.push_back(-9999);
      m_y_flt.push_back(-9999);
      m_z_flt.push_back(-9999);
      m_py_flt.push_back(-9999);
      m_pz_flt.push_back(-9999);
      m_pT_flt.push_back(-9999);
      m_eta_flt.push_back(-9999);
      m_chi2.push_back(-9999);
    }
  
    bool smoothed = false;
    if (state.hasSmoothed())
    {
      smoothed = true;
      m_nSmoothed++;

      auto parameter = state.smoothed();
      auto covariance = state.smoothedCovariance();

      m_eLOC0_smt.push_back(parameter[Acts::eBoundLoc0]);
      m_eLOC1_smt.push_back(parameter[Acts::eBoundLoc1]);
      m_ePHI_smt.push_back(parameter[Acts::eBoundPhi]);
      m_eTHETA_smt.push_back(parameter[Acts::eBoundTheta]);
      m_eQOP_smt.push_back(parameter[Acts::eBoundQOverP]);
      m_eT_smt.push_back(parameter[Acts::eBoundTime]);

      m_res_eLOC0_smt.push_back(parameter[Acts::eBoundLoc0] -
                                truthLOC0);
      m_res_eLOC1_smt.push_back(parameter[Acts::eBoundLoc1] -
                                truthLOC1);
      m_res_ePHI_smt.push_back(parameter[Acts::eBoundPhi] -
                               truthPHI);
      m_res_eTHETA_smt.push_back(parameter[Acts::eBoundTheta] -
                                 truthTHETA);
      m_res_eQOP_smt.push_back(parameter[Acts::eBoundQOverP] -
                               truthQOP);
      m_res_eT_smt.push_back(parameter[Acts::eBoundTime] -
                             truthTIME);

      /// Smoothed parameter uncertainties
      m_err_eLOC0_smt.push_back(
          sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_err_eLOC1_smt.push_back(
          sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_err_ePHI_smt.push_back(
          sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
      m_err_eTHETA_smt.push_back(
          sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
      m_err_eQOP_smt.push_back(
          sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
      m_err_eT_smt.push_back(
          sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

      /// Smoothed parameter pulls
      m_pull_eLOC0_smt.push_back(
          (parameter[Acts::eBoundLoc0] - truthLOC0) /
          sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0)));
      m_pull_eLOC1_smt.push_back(
          (parameter[Acts::eBoundLoc1] - truthLOC1) /
          sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1)));
      m_pull_ePHI_smt.push_back(
          (parameter[Acts::eBoundPhi] - truthPHI) /
          sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi)));
      m_pull_eTHETA_smt.push_back(
          (parameter[Acts::eBoundTheta] - truthTHETA) /
          sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta)));
      m_pull_eQOP_smt.push_back(
          (parameter[Acts::eBoundQOverP] - truthQOP) /
          sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP)));
      m_pull_eT_smt.push_back(
          (parameter[Acts::eBoundTime] - truthTIME) /
          sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime)));

      Acts::FreeVector freeparams = Acts::detail::transformBoundToFreeParameters(surface, m_tGeometry->geometry().getGeoContext(), parameter);
      
      /// Other smoothed parameter info
      m_x_smt.push_back(freeparams[Acts::eFreePos0]);
      m_y_smt.push_back(freeparams[Acts::eFreePos1]);
      m_z_smt.push_back(freeparams[Acts::eFreePos2]);
      auto p = std::abs(1/freeparams[Acts::eFreeQOverP]);
      m_px_smt.push_back(p * freeparams[Acts::eFreeDir0]);
      m_py_smt.push_back(p * freeparams[Acts::eFreeDir1]);
      m_pz_smt.push_back(p * freeparams[Acts::eFreeDir2]);
      m_pT_smt.push_back(sqrt(p * std::hypot(freeparams[Acts::eFreeDir0],
					     freeparams[Acts::eFreeDir1])));
      m_eta_smt.push_back(eta(freeparams.segment<3>(Acts::eFreeDir0)));

    }
    else
    {
      /// Push bad values if no smoothed parameter
      m_eLOC0_smt.push_back(-9999);
      m_eLOC1_smt.push_back(-9999);
      m_ePHI_smt.push_back(-9999);
      m_eTHETA_smt.push_back(-9999);
      m_eQOP_smt.push_back(-9999);
      m_eT_smt.push_back(-9999);
      m_res_eLOC0_smt.push_back(-9999);
      m_res_eLOC1_smt.push_back(-9999);
      m_res_ePHI_smt.push_back(-9999);
      m_res_eTHETA_smt.push_back(-9999);
      m_res_eQOP_smt.push_back(-9999);
      m_res_eT_smt.push_back(-9999);
      m_err_eLOC0_smt.push_back(-9999);
      m_err_eLOC1_smt.push_back(-9999);
      m_err_ePHI_smt.push_back(-9999);
      m_err_eTHETA_smt.push_back(-9999);
      m_err_eQOP_smt.push_back(-9999);
      m_err_eT_smt.push_back(-9999);
      m_pull_eLOC0_smt.push_back(-9999);
      m_pull_eLOC1_smt.push_back(-9999);
      m_pull_ePHI_smt.push_back(-9999);
      m_pull_eTHETA_smt.push_back(-9999);
      m_pull_eQOP_smt.push_back(-9999);
      m_pull_eT_smt.push_back(-9999);
      m_x_smt.push_back(-9999);
      m_y_smt.push_back(-9999);
      m_z_smt.push_back(-9999);
      m_px_smt.push_back(-9999);
      m_py_smt.push_back(-9999);
      m_pz_smt.push_back(-9999);
      m_pT_smt.push_back(-9999);
      m_eta_smt.push_back(-9999);
    }

    /// Save whether or not states had various KF steps
    m_prt.push_back(predicted);
    m_flt.push_back(filtered);
    m_smt.push_back(smoothed);

    return true; }  /// Finish lambda function
  );                       /// Finish multi trajectory visitBackwards call

  if (m_verbosity > 2)
    std::cout << "Finished track states" << std::endl;

  return;
}

Acts::Vector3 ActsEvaluator::getGlobalTruthHit(TrkrDefs::cluskey cluskey,
                                               float& _gt)
{
  SvtxClusterEval* clustereval = m_svtxEvalStack->get_cluster_eval();

  const auto [truth_ckey, truth_cluster] = clustereval->max_truth_cluster_by_energy(cluskey);

  float gx = -9999;
  float gy = -9999;
  float gz = -9999;
  float gt = -9999;

  if (truth_cluster)
  {
    gx = truth_cluster->getX();
    gy = truth_cluster->getY();
    gz = truth_cluster->getZ();
  }

  /// Convert to acts units of mm
  gx *= Acts::UnitConstants::cm;
  gy *= Acts::UnitConstants::cm;
  gz *= Acts::UnitConstants::cm;

  Acts::Vector3 globalPos(gx, gy, gz);
  _gt = gt;
  return globalPos;
}

//___________________________________________________________________________________
Surface ActsEvaluator::getSurface(TrkrDefs::cluskey cluskey, TrkrCluster* cluster)
{
  return m_tGeometry->maps().getSurface(cluskey, cluster);
}

void ActsEvaluator::fillProtoTrack(const TrackSeed* seed)
{
  if (m_verbosity > 2)
  {
    std::cout << "Filling proto track seed quantities" << std::endl;
  }

  unsigned int tpcid = seed->get_tpc_seed_index();
  unsigned int siid = seed->get_silicon_seed_index();

  auto siseed = m_siliconSeeds->get(siid);
  auto tpcseed = m_tpcSeeds->get(tpcid);
  if(!tpcseed) return;

  Acts::Vector3 position = Acts::Vector3::Zero();

  if (siseed)
  {
    position(0) = siseed->get_x() * Acts::UnitConstants::cm;
    position(1) = siseed->get_y() * Acts::UnitConstants::cm;
    position(2) = siseed->get_z() * Acts::UnitConstants::cm;
  }
  else
  {
    position(0) = tpcseed->get_x() * Acts::UnitConstants::cm;
    position(1) = tpcseed->get_y() * Acts::UnitConstants::cm;
    position(2) = tpcseed->get_z() * Acts::UnitConstants::cm;
  }

  Acts::Vector3 momentum(tpcseed->get_px(m_clusterContainer,
                                         m_tGeometry),
                         tpcseed->get_py(m_clusterContainer,
                                         m_tGeometry),
                         tpcseed->get_pz());

  m_protoTrackPx = momentum(0);
  m_protoTrackPy = momentum(1);
  m_protoTrackPz = momentum(2);
  m_protoTrackX = position(0);
  m_protoTrackY = position(1);
  m_protoTrackZ = position(2);

  for (auto svtxseed : {siseed, tpcseed})
  {
    //! protect against tpc only tracks
    if (!svtxseed)
    {
      continue;
    }
    for (auto clusIter = svtxseed->begin_cluster_keys();
         clusIter != svtxseed->end_cluster_keys();
         ++clusIter)
    {
      auto key = *clusIter;
      auto cluster = m_clusterContainer->findCluster(key);

      /// Get source link global position
      Acts::Vector2 loc(cluster->getLocalX(),
                        cluster->getLocalY());

      if (TrkrDefs::getTrkrId(key) == TrkrDefs::TrkrId::tpcId)
      {
        // must convert local Y from cluster average time of arival to local cluster z position
        double drift_velocity = m_tGeometry->get_drift_velocity();
        double zdriftlength = cluster->getLocalY() * drift_velocity;
        double surfCenterZ = 52.89;                // 52.89 is where G4 thinks the surface center is
        double zloc = surfCenterZ - zdriftlength;  // converts z drift length to local z position in the TPC in north
        unsigned int side = TpcDefs::getSide(key);
        if (side == 0) zloc = -zloc;
        loc(1) = zloc * Acts::UnitConstants::cm;
      }

      Acts::Vector3 mom(0, 0, 0);
      Acts::Vector3 globalPos = m_tGeometry->getGlobalPosition(key, cluster) * Acts::UnitConstants::cm;

      m_SLx.push_back(globalPos(0));
      m_SLy.push_back(globalPos(1));
      m_SLz.push_back(globalPos(2));
      m_SL_lx.push_back(loc(0));
      m_SL_ly.push_back(loc(1));

      /// Get corresponding truth hit position
      float gt = -9999;

      Acts::Vector3 globalTruthPos = getGlobalTruthHit(key, gt);

      float gx = globalTruthPos(0);
      float gy = globalTruthPos(1);
      float gz = globalTruthPos(2);

      /// Get local truth position
      const float r = sqrt(gx * gx + gy * gy + gz * gz);
      Acts::Vector3 globalTruthUnitDir(gx / r, gy / r, gz / r);

      auto surf = getSurface(key, cluster);

      auto truthLocal = (*surf).globalToLocal(m_tGeometry->geometry().getGeoContext(),
                                              globalTruthPos,
                                              globalTruthUnitDir);

      if (truthLocal.ok())
      {
        Acts::Vector2 truthLocalVec = truthLocal.value();

        m_t_SL_lx.push_back(truthLocalVec(0));
        m_t_SL_ly.push_back(truthLocalVec(1));
        m_t_SL_gx.push_back(gx);
        m_t_SL_gy.push_back(gy);
        m_t_SL_gz.push_back(gz);
      }
      else
      {
        Acts::Vector3 loct = (*surf).transform(m_tGeometry->geometry().getGeoContext()).inverse() * globalTruthPos;

        m_t_SL_lx.push_back(loct(0));
        m_t_SL_ly.push_back(loct(1));
        m_t_SL_gx.push_back(gx);
        m_t_SL_gy.push_back(gy);
        m_t_SL_gz.push_back(gz);
      }
    }
  }
  if (m_verbosity > 2)
  {
    std::cout << "Filled proto track" << std::endl;
  }
}

void ActsEvaluator::fillFittedTrackParams(const Trajectory::IndexedParameters& paramsMap,
                                          const size_t& trackTip)
{
  m_hasFittedParams = false;

  if (m_verbosity > 2)
  {
    std::cout << "Filling fitted track parameters" << std::endl;
  }

  /// If it has track parameters, fill the values
  if (paramsMap.find(trackTip) != paramsMap.end())
  {
    m_hasFittedParams = true;
    const auto& boundParam = paramsMap.find(trackTip)->second;
    const auto& parameter = boundParam.parameters();
    const auto& covariance = *boundParam.covariance();
    m_charge_fit = boundParam.charge();
    m_eLOC0_fit = parameter[Acts::eBoundLoc0];
    m_eLOC1_fit = parameter[Acts::eBoundLoc1];
    m_ePHI_fit = parameter[Acts::eBoundPhi];
    m_eTHETA_fit = parameter[Acts::eBoundTheta];
    m_eQOP_fit = parameter[Acts::eBoundQOverP];
    m_eT_fit = parameter[Acts::eBoundTime];
    m_err_eLOC0_fit =
        sqrt(covariance(Acts::eBoundLoc0, Acts::eBoundLoc0));
    m_err_eLOC1_fit =
        sqrt(covariance(Acts::eBoundLoc1, Acts::eBoundLoc1));
    m_err_ePHI_fit = sqrt(covariance(Acts::eBoundPhi, Acts::eBoundPhi));
    m_err_eTHETA_fit =
        sqrt(covariance(Acts::eBoundTheta, Acts::eBoundTheta));
    m_err_eQOP_fit = sqrt(covariance(Acts::eBoundQOverP, Acts::eBoundQOverP));
    m_err_eT_fit = sqrt(covariance(Acts::eBoundTime, Acts::eBoundTime));

    m_px_fit = boundParam.momentum()(0);
    m_py_fit = boundParam.momentum()(1);
    m_pz_fit = boundParam.momentum()(2);
    m_x_fit = boundParam.position(m_tGeometry->geometry().getGeoContext())(0);
    m_y_fit = boundParam.position(m_tGeometry->geometry().getGeoContext())(1);
    m_z_fit = boundParam.position(m_tGeometry->geometry().getGeoContext())(2);

    return;
  }

  /// Otherwise mark it as a bad fit
  m_eLOC0_fit = -9999;
  m_eLOC1_fit = -9999;
  m_ePHI_fit = -9999;
  m_eTHETA_fit = -9999;
  m_eQOP_fit = -9999;
  m_eT_fit = -9999;
  m_charge_fit = -9999;
  m_err_eLOC0_fit = -9999;
  m_err_eLOC1_fit = -9999;
  m_err_ePHI_fit = -9999;
  m_err_eTHETA_fit = -9999;
  m_err_eQOP_fit = -9999;
  m_err_eT_fit = -9999;
  m_px_fit = -9999;
  m_py_fit = -9999;
  m_pz_fit = -9999;
  m_x_fit = -9999;
  m_y_fit = -9999;
  m_z_fit = -9999;

  if (m_verbosity > 2)
  {
    std::cout << "Finished fitted track params" << std::endl;
  }
  return;
}

void ActsEvaluator::fillG4Particle(PHG4Particle* part)
{
  SvtxTruthEval* trutheval = m_svtxEvalStack->get_truth_eval();

  if (part)
  {
    m_t_barcode = part->get_track_id();
    const auto pid = part->get_pid();
    m_t_charge = pid < 0 ? -1 : 1;
    const auto vtx = trutheval->get_vertex(part);
    m_t_vx = vtx->get_x() * Acts::UnitConstants::cm;
    m_t_vy = vtx->get_y() * Acts::UnitConstants::cm;
    m_t_vz = vtx->get_z() * Acts::UnitConstants::cm;
    if (m_verbosity > 1)
      std::cout << "truth vertex : (" << m_t_vx << ", " << m_t_vy
                << ", " << m_t_vz << ")" << std::endl;
    m_t_px = part->get_px();
    m_t_py = part->get_py();
    m_t_pz = part->get_pz();
    const double p = sqrt(m_t_px * m_t_px + m_t_py * m_t_py + m_t_pz * m_t_pz);
    m_t_theta = acos(m_t_pz / p);
    m_t_phi = atan(m_t_py / m_t_px);
    m_t_pT = sqrt(m_t_px * m_t_px + m_t_py * m_t_py);
    m_t_eta = atanh(m_t_pz / p);

    return;
  }

  /// If particle doesn't exist, just fill with -999
  m_t_barcode = -9999;
  m_t_charge = -9999;
  m_t_vx = -9999;
  m_t_vy = -9999;
  m_t_vz = -9999;
  m_t_px = -9999;
  m_t_py = -9999;
  m_t_pz = -9999;
  m_t_theta = -9999;
  m_t_phi = -9999;
  m_t_pT = -9999;
  m_t_eta = -9999;

  return;
}

int ActsEvaluator::getNodes(PHCompositeNode* topNode)
{
  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");

  if (!m_tpcSeeds or !m_siliconSeeds)
  {
    std::cout << PHWHERE << "Seed containers not found, cannot continue!"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!m_truthInfo)
  {
    std::cout << PHWHERE << "PHG4TruthInfoContainer not found, cannot continue!"
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No Acts Tracking geometry on node tree. Bailing"
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!m_trackMap)
  {
    std::cout << PHWHERE << "No SvtxTrackMap on node tree. Bailing."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_actsProtoTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SeedTrackMap");
  if (!m_actsProtoTrackMap)
  {
    std::cout << PHWHERE << "No Acts proto tracks on node tree. Bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE << "No clusters, bailing"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void ActsEvaluator::clearTrackVariables()
{
  m_t_x.clear();
  m_t_y.clear();
  m_t_z.clear();
  m_t_r.clear();
  m_t_dx.clear();
  m_t_dy.clear();
  m_t_dz.clear();
  m_t_eLOC0.clear();
  m_t_eLOC1.clear();
  m_t_ePHI.clear();
  m_t_eTHETA.clear();
  m_t_eQOP.clear();
  m_t_eT.clear();

  m_volumeID.clear();
  m_layerID.clear();
  m_moduleID.clear();
  m_lx_hit.clear();
  m_ly_hit.clear();
  m_x_hit.clear();
  m_y_hit.clear();
  m_z_hit.clear();
  m_res_x_hit.clear();
  m_res_y_hit.clear();
  m_err_x_hit.clear();
  m_err_y_hit.clear();
  m_pull_x_hit.clear();
  m_pull_y_hit.clear();
  m_dim_hit.clear();

  m_prt.clear();
  m_eLOC0_prt.clear();
  m_eLOC1_prt.clear();
  m_ePHI_prt.clear();
  m_eTHETA_prt.clear();
  m_eQOP_prt.clear();
  m_eT_prt.clear();
  m_res_eLOC0_prt.clear();
  m_res_eLOC1_prt.clear();
  m_res_ePHI_prt.clear();
  m_res_eTHETA_prt.clear();
  m_res_eQOP_prt.clear();
  m_res_eT_prt.clear();
  m_err_eLOC0_prt.clear();
  m_err_eLOC1_prt.clear();
  m_err_ePHI_prt.clear();
  m_err_eTHETA_prt.clear();
  m_err_eQOP_prt.clear();
  m_err_eT_prt.clear();
  m_pull_eLOC0_prt.clear();
  m_pull_eLOC1_prt.clear();
  m_pull_ePHI_prt.clear();
  m_pull_eTHETA_prt.clear();
  m_pull_eQOP_prt.clear();
  m_pull_eT_prt.clear();
  m_x_prt.clear();
  m_y_prt.clear();
  m_z_prt.clear();
  m_px_prt.clear();
  m_py_prt.clear();
  m_pz_prt.clear();
  m_eta_prt.clear();
  m_pT_prt.clear();

  m_flt.clear();
  m_eLOC0_flt.clear();
  m_eLOC1_flt.clear();
  m_ePHI_flt.clear();
  m_eTHETA_flt.clear();
  m_eQOP_flt.clear();
  m_eT_flt.clear();
  m_res_eLOC0_flt.clear();
  m_res_eLOC1_flt.clear();
  m_res_ePHI_flt.clear();
  m_res_eTHETA_flt.clear();
  m_res_eQOP_flt.clear();
  m_res_eT_flt.clear();
  m_err_eLOC0_flt.clear();
  m_err_eLOC1_flt.clear();
  m_err_ePHI_flt.clear();
  m_err_eTHETA_flt.clear();
  m_err_eQOP_flt.clear();
  m_err_eT_flt.clear();
  m_pull_eLOC0_flt.clear();
  m_pull_eLOC1_flt.clear();
  m_pull_ePHI_flt.clear();
  m_pull_eTHETA_flt.clear();
  m_pull_eQOP_flt.clear();
  m_pull_eT_flt.clear();
  m_x_flt.clear();
  m_y_flt.clear();
  m_z_flt.clear();
  m_px_flt.clear();
  m_py_flt.clear();
  m_pz_flt.clear();
  m_eta_flt.clear();
  m_pT_flt.clear();
  m_chi2.clear();

  m_smt.clear();
  m_eLOC0_smt.clear();
  m_eLOC1_smt.clear();
  m_ePHI_smt.clear();
  m_eTHETA_smt.clear();
  m_eQOP_smt.clear();
  m_eT_smt.clear();
  m_res_eLOC0_smt.clear();
  m_res_eLOC1_smt.clear();
  m_res_ePHI_smt.clear();
  m_res_eTHETA_smt.clear();
  m_res_eQOP_smt.clear();
  m_res_eT_smt.clear();
  m_err_eLOC0_smt.clear();
  m_err_eLOC1_smt.clear();
  m_err_ePHI_smt.clear();
  m_err_eTHETA_smt.clear();
  m_err_eQOP_smt.clear();
  m_err_eT_smt.clear();
  m_pull_eLOC0_smt.clear();
  m_pull_eLOC1_smt.clear();
  m_pull_ePHI_smt.clear();
  m_pull_eTHETA_smt.clear();
  m_pull_eQOP_smt.clear();
  m_pull_eT_smt.clear();
  m_x_smt.clear();
  m_y_smt.clear();
  m_z_smt.clear();
  m_px_smt.clear();
  m_py_smt.clear();
  m_pz_smt.clear();
  m_eta_smt.clear();
  m_pT_smt.clear();

  m_SLx.clear();
  m_SLy.clear();
  m_SLz.clear();
  m_SL_lx.clear();
  m_SL_ly.clear();
  m_t_SL_lx.clear();
  m_t_SL_ly.clear();
  m_t_SL_gx.clear();
  m_t_SL_gy.clear();
  m_t_SL_gz.clear();

  m_protoTrackPx = -9999.;
  m_protoTrackPy = -9999.;
  m_protoTrackPz = -9999.;
  m_protoTrackX = -9999.;
  m_protoTrackY = -9999.;
  m_protoTrackZ = -9999.;
  m_protoD0Cov = -9999.;
  m_protoZ0Cov = -9999.;
  m_protoPhiCov = -9999.;
  m_protoThetaCov = -9999.;
  m_protoQopCov = -9999.;

  return;
}

void ActsEvaluator::initializeTree()
{
  m_trackFile = new TFile(m_filename.c_str(), "RECREATE");

  m_trackTree = new TTree("tracktree", "A tree with Acts KF track information");

  m_trackTree->Branch("event_nr", &m_eventNr);
  m_trackTree->Branch("traj_nr", &m_trajNr);
  m_trackTree->Branch("track_nr", &m_trackNr);
  m_trackTree->Branch("t_barcode", &m_t_barcode, "t_barcode/l");
  m_trackTree->Branch("t_charge", &m_t_charge);
  m_trackTree->Branch("t_time", &m_t_time);
  m_trackTree->Branch("t_vx", &m_t_vx);
  m_trackTree->Branch("t_vy", &m_t_vy);
  m_trackTree->Branch("t_vz", &m_t_vz);
  m_trackTree->Branch("t_px", &m_t_px);
  m_trackTree->Branch("t_py", &m_t_py);
  m_trackTree->Branch("t_pz", &m_t_pz);
  m_trackTree->Branch("t_theta", &m_t_theta);
  m_trackTree->Branch("t_phi", &m_t_phi);
  m_trackTree->Branch("t_eta", &m_t_eta);
  m_trackTree->Branch("t_pT", &m_t_pT);

  m_trackTree->Branch("t_x", &m_t_x);
  m_trackTree->Branch("t_y", &m_t_y);
  m_trackTree->Branch("t_z", &m_t_z);
  m_trackTree->Branch("t_r", &m_t_r);
  m_trackTree->Branch("t_dx", &m_t_dx);
  m_trackTree->Branch("t_dy", &m_t_dy);
  m_trackTree->Branch("t_dz", &m_t_dz);
  m_trackTree->Branch("t_eLOC0", &m_t_eLOC0);
  m_trackTree->Branch("t_eLOC1", &m_t_eLOC1);
  m_trackTree->Branch("t_ePHI", &m_t_ePHI);
  m_trackTree->Branch("t_eTHETA", &m_t_eTHETA);
  m_trackTree->Branch("t_eQOP", &m_t_eQOP);
  m_trackTree->Branch("t_eT", &m_t_eT);

  m_trackTree->Branch("hasFittedParams", &m_hasFittedParams);
  m_trackTree->Branch("chi2_fit", &m_chi2_fit);
  m_trackTree->Branch("quality", &m_quality);
  m_trackTree->Branch("ndf_fit", &m_ndf_fit);
  m_trackTree->Branch("eLOC0_fit", &m_eLOC0_fit);
  m_trackTree->Branch("eLOC1_fit", &m_eLOC1_fit);
  m_trackTree->Branch("ePHI_fit", &m_ePHI_fit);
  m_trackTree->Branch("eTHETA_fit", &m_eTHETA_fit);
  m_trackTree->Branch("eQOP_fit", &m_eQOP_fit);
  m_trackTree->Branch("eT_fit", &m_eT_fit);
  m_trackTree->Branch("charge_fit", &m_charge_fit);
  m_trackTree->Branch("err_eLOC0_fit", &m_err_eLOC0_fit);
  m_trackTree->Branch("err_eLOC1_fit", &m_err_eLOC1_fit);
  m_trackTree->Branch("err_ePHI_fit", &m_err_ePHI_fit);
  m_trackTree->Branch("err_eTHETA_fit", &m_err_eTHETA_fit);
  m_trackTree->Branch("err_eQOP_fit", &m_err_eQOP_fit);
  m_trackTree->Branch("err_eT_fit", &m_err_eT_fit);
  m_trackTree->Branch("g_px_fit", &m_px_fit);
  m_trackTree->Branch("g_py_fit", &m_py_fit);
  m_trackTree->Branch("g_pz_fit", &m_pz_fit);
  m_trackTree->Branch("g_x_fit", &m_x_fit);
  m_trackTree->Branch("g_y_fit", &m_y_fit);
  m_trackTree->Branch("g_z_fit", &m_z_fit);
  m_trackTree->Branch("g_dca3Dxy_fit", &m_dca3Dxy);
  m_trackTree->Branch("g_dca3Dz_fit", &m_dca3Dz);
  m_trackTree->Branch("g_dca3Dxy_cov", &m_dca3DxyCov);
  m_trackTree->Branch("g_dca3Dz_cov", &m_dca3DzCov);

  m_trackTree->Branch("g_protoTrackPx", &m_protoTrackPx);
  m_trackTree->Branch("g_protoTrackPy", &m_protoTrackPy);
  m_trackTree->Branch("g_protoTrackPz", &m_protoTrackPz);
  m_trackTree->Branch("g_protoTrackX", &m_protoTrackX);
  m_trackTree->Branch("g_protoTrackY", &m_protoTrackY);
  m_trackTree->Branch("g_protoTrackZ", &m_protoTrackZ);
  m_trackTree->Branch("g_protoTrackD0Cov", &m_protoD0Cov);
  m_trackTree->Branch("g_protoTrackZ0Cov", &m_protoZ0Cov);
  m_trackTree->Branch("g_protoTrackPhiCov", &m_protoPhiCov);
  m_trackTree->Branch("g_protoTrackThetaCov", &m_protoThetaCov);
  m_trackTree->Branch("g_protoTrackQopCov", &m_protoQopCov);
  m_trackTree->Branch("g_SLx", &m_SLx);
  m_trackTree->Branch("g_SLy", &m_SLy);
  m_trackTree->Branch("g_SLz", &m_SLz);
  m_trackTree->Branch("l_SLx", &m_SL_lx);
  m_trackTree->Branch("l_SLy", &m_SL_ly);
  m_trackTree->Branch("t_SL_lx", &m_t_SL_lx);
  m_trackTree->Branch("t_SL_ly", &m_t_SL_ly);
  m_trackTree->Branch("t_SL_gx", &m_t_SL_gx);
  m_trackTree->Branch("t_SL_gy", &m_t_SL_gy);
  m_trackTree->Branch("t_SL_gz", &m_t_SL_gz);

  m_trackTree->Branch("nSharedHits", &m_nSharedHits);
  m_trackTree->Branch("nHoles", &m_nHoles);
  m_trackTree->Branch("nOutliers", &m_nOutliers);
  m_trackTree->Branch("nStates", &m_nStates);
  m_trackTree->Branch("nMeasurements", &m_nMeasurements);
  m_trackTree->Branch("volume_id", &m_volumeID);
  m_trackTree->Branch("layer_id", &m_layerID);
  m_trackTree->Branch("module_id", &m_moduleID);
  m_trackTree->Branch("l_x_hit", &m_lx_hit);
  m_trackTree->Branch("l_y_hit", &m_ly_hit);
  m_trackTree->Branch("g_x_hit", &m_x_hit);
  m_trackTree->Branch("g_y_hit", &m_y_hit);
  m_trackTree->Branch("g_z_hit", &m_z_hit);
  m_trackTree->Branch("res_x_hit", &m_res_x_hit);
  m_trackTree->Branch("res_y_hit", &m_res_y_hit);
  m_trackTree->Branch("err_x_hit", &m_err_x_hit);
  m_trackTree->Branch("err_y_hit", &m_err_y_hit);
  m_trackTree->Branch("pull_x_hit", &m_pull_x_hit);
  m_trackTree->Branch("pull_y_hit", &m_pull_y_hit);
  m_trackTree->Branch("dim_hit", &m_dim_hit);

  m_trackTree->Branch("nPredicted", &m_nPredicted);
  m_trackTree->Branch("predicted", &m_prt);
  m_trackTree->Branch("eLOC0_prt", &m_eLOC0_prt);
  m_trackTree->Branch("eLOC1_prt", &m_eLOC1_prt);
  m_trackTree->Branch("ePHI_prt", &m_ePHI_prt);
  m_trackTree->Branch("eTHETA_prt", &m_eTHETA_prt);
  m_trackTree->Branch("eQOP_prt", &m_eQOP_prt);
  m_trackTree->Branch("eT_prt", &m_eT_prt);
  m_trackTree->Branch("res_eLOC0_prt", &m_res_eLOC0_prt);
  m_trackTree->Branch("res_eLOC1_prt", &m_res_eLOC1_prt);
  m_trackTree->Branch("res_ePHI_prt", &m_res_ePHI_prt);
  m_trackTree->Branch("res_eTHETA_prt", &m_res_eTHETA_prt);
  m_trackTree->Branch("res_eQOP_prt", &m_res_eQOP_prt);
  m_trackTree->Branch("res_eT_prt", &m_res_eT_prt);
  m_trackTree->Branch("err_eLOC0_prt", &m_err_eLOC0_prt);
  m_trackTree->Branch("err_eLOC1_prt", &m_err_eLOC1_prt);
  m_trackTree->Branch("err_ePHI_prt", &m_err_ePHI_prt);
  m_trackTree->Branch("err_eTHETA_prt", &m_err_eTHETA_prt);
  m_trackTree->Branch("err_eQOP_prt", &m_err_eQOP_prt);
  m_trackTree->Branch("err_eT_prt", &m_err_eT_prt);
  m_trackTree->Branch("pull_eLOC0_prt", &m_pull_eLOC0_prt);
  m_trackTree->Branch("pull_eLOC1_prt", &m_pull_eLOC1_prt);
  m_trackTree->Branch("pull_ePHI_prt", &m_pull_ePHI_prt);
  m_trackTree->Branch("pull_eTHETA_prt", &m_pull_eTHETA_prt);
  m_trackTree->Branch("pull_eQOP_prt", &m_pull_eQOP_prt);
  m_trackTree->Branch("pull_eT_prt", &m_pull_eT_prt);
  m_trackTree->Branch("g_x_prt", &m_x_prt);
  m_trackTree->Branch("g_y_prt", &m_y_prt);
  m_trackTree->Branch("g_z_prt", &m_z_prt);
  m_trackTree->Branch("px_prt", &m_px_prt);
  m_trackTree->Branch("py_prt", &m_py_prt);
  m_trackTree->Branch("pz_prt", &m_pz_prt);
  m_trackTree->Branch("eta_prt", &m_eta_prt);
  m_trackTree->Branch("pT_prt", &m_pT_prt);

  m_trackTree->Branch("nFiltered", &m_nFiltered);
  m_trackTree->Branch("filtered", &m_flt);
  m_trackTree->Branch("eLOC0_flt", &m_eLOC0_flt);
  m_trackTree->Branch("eLOC1_flt", &m_eLOC1_flt);
  m_trackTree->Branch("ePHI_flt", &m_ePHI_flt);
  m_trackTree->Branch("eTHETA_flt", &m_eTHETA_flt);
  m_trackTree->Branch("eQOP_flt", &m_eQOP_flt);
  m_trackTree->Branch("eT_flt", &m_eT_flt);
  m_trackTree->Branch("res_eLOC0_flt", &m_res_eLOC0_flt);
  m_trackTree->Branch("res_eLOC1_flt", &m_res_eLOC1_flt);
  m_trackTree->Branch("res_ePHI_flt", &m_res_ePHI_flt);
  m_trackTree->Branch("res_eTHETA_flt", &m_res_eTHETA_flt);
  m_trackTree->Branch("res_eQOP_flt", &m_res_eQOP_flt);
  m_trackTree->Branch("res_eT_flt", &m_res_eT_flt);
  m_trackTree->Branch("err_eLOC0_flt", &m_err_eLOC0_flt);
  m_trackTree->Branch("err_eLOC1_flt", &m_err_eLOC1_flt);
  m_trackTree->Branch("err_ePHI_flt", &m_err_ePHI_flt);
  m_trackTree->Branch("err_eTHETA_flt", &m_err_eTHETA_flt);
  m_trackTree->Branch("err_eQOP_flt", &m_err_eQOP_flt);
  m_trackTree->Branch("err_eT_flt", &m_err_eT_flt);
  m_trackTree->Branch("pull_eLOC0_flt", &m_pull_eLOC0_flt);
  m_trackTree->Branch("pull_eLOC1_flt", &m_pull_eLOC1_flt);
  m_trackTree->Branch("pull_ePHI_flt", &m_pull_ePHI_flt);
  m_trackTree->Branch("pull_eTHETA_flt", &m_pull_eTHETA_flt);
  m_trackTree->Branch("pull_eQOP_flt", &m_pull_eQOP_flt);
  m_trackTree->Branch("pull_eT_flt", &m_pull_eT_flt);
  m_trackTree->Branch("g_x_flt", &m_x_flt);
  m_trackTree->Branch("g_y_flt", &m_y_flt);
  m_trackTree->Branch("g_z_flt", &m_z_flt);
  m_trackTree->Branch("px_flt", &m_px_flt);
  m_trackTree->Branch("py_flt", &m_py_flt);
  m_trackTree->Branch("pz_flt", &m_pz_flt);
  m_trackTree->Branch("eta_flt", &m_eta_flt);
  m_trackTree->Branch("pT_flt", &m_pT_flt);
  m_trackTree->Branch("chi2", &m_chi2);

  m_trackTree->Branch("nSmoothed", &m_nSmoothed);
  m_trackTree->Branch("smoothed", &m_smt);
  m_trackTree->Branch("eLOC0_smt", &m_eLOC0_smt);
  m_trackTree->Branch("eLOC1_smt", &m_eLOC1_smt);
  m_trackTree->Branch("ePHI_smt", &m_ePHI_smt);
  m_trackTree->Branch("eTHETA_smt", &m_eTHETA_smt);
  m_trackTree->Branch("eQOP_smt", &m_eQOP_smt);
  m_trackTree->Branch("eT_smt", &m_eT_smt);
  m_trackTree->Branch("res_eLOC0_smt", &m_res_eLOC0_smt);
  m_trackTree->Branch("res_eLOC1_smt", &m_res_eLOC1_smt);
  m_trackTree->Branch("res_ePHI_smt", &m_res_ePHI_smt);
  m_trackTree->Branch("res_eTHETA_smt", &m_res_eTHETA_smt);
  m_trackTree->Branch("res_eQOP_smt", &m_res_eQOP_smt);
  m_trackTree->Branch("res_eT_smt", &m_res_eT_smt);
  m_trackTree->Branch("err_eLOC0_smt", &m_err_eLOC0_smt);
  m_trackTree->Branch("err_eLOC1_smt", &m_err_eLOC1_smt);
  m_trackTree->Branch("err_ePHI_smt", &m_err_ePHI_smt);
  m_trackTree->Branch("err_eTHETA_smt", &m_err_eTHETA_smt);
  m_trackTree->Branch("err_eQOP_smt", &m_err_eQOP_smt);
  m_trackTree->Branch("err_eT_smt", &m_err_eT_smt);
  m_trackTree->Branch("pull_eLOC0_smt", &m_pull_eLOC0_smt);
  m_trackTree->Branch("pull_eLOC1_smt", &m_pull_eLOC1_smt);
  m_trackTree->Branch("pull_ePHI_smt", &m_pull_ePHI_smt);
  m_trackTree->Branch("pull_eTHETA_smt", &m_pull_eTHETA_smt);
  m_trackTree->Branch("pull_eQOP_smt", &m_pull_eQOP_smt);
  m_trackTree->Branch("pull_eT_smt", &m_pull_eT_smt);
  m_trackTree->Branch("g_x_smt", &m_x_smt);
  m_trackTree->Branch("g_y_smt", &m_y_smt);
  m_trackTree->Branch("g_z_smt", &m_z_smt);
  m_trackTree->Branch("px_smt", &m_px_smt);
  m_trackTree->Branch("py_smt", &m_py_smt);
  m_trackTree->Branch("pz_smt", &m_pz_smt);
  m_trackTree->Branch("eta_smt", &m_eta_smt);
  m_trackTree->Branch("pT_smt", &m_pT_smt);
}
