#include "ActsEvaluator.h"

/// General fun4all and subsysreco includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

/// Tracking includes
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxEvaluator.h>
#include <g4eval/SvtxTrackEval.h>

#include <g4main/PHG4VtxPoint.h>

#include <Acts/EventData/MultiTrajectoryHelpers.hpp>
#include <Acts/Utilities/Units.hpp>

#include <TFile.h>
#include <TTree.h>

ActsEvaluator::ActsEvaluator(const std::string &name,
                                       SvtxEvaluator *svtxEvaluator)
  : SubsysReco(name)
  , m_svtxEvaluator(svtxEvaluator)
  , m_truthInfo(nullptr)
  , m_trackMap(nullptr)
  , m_svtxEvalStack(nullptr)
  , m_actsFitResults(nullptr)
  , m_hitIdClusKey(nullptr)
  , m_tGeometry(nullptr)
{
}

ActsEvaluator::~ActsEvaluator()
{
}

int ActsEvaluator::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Starting ActsEvaluator::Init" << std::endl;
  }

  initializeTree();

  if (Verbosity() > 1)
  {
    std::cout << "Finished ActsEvaluator::Init" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int ActsEvaluator::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << "Starting ActsEvaluator at event " << m_eventNr
              << std::endl;
  }
  
  if (getNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
    return Fun4AllReturnCodes::ABORTEVENT;

  m_svtxEvalStack = new SvtxEvalStack(topNode);

  m_svtxEvalStack->next_event(topNode);

  SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();

  std::map<const unsigned int, Trajectory>::iterator trackIter;
  int iTrack = 0;

  for (trackIter = m_actsFitResults->begin();
       trackIter != m_actsFitResults->end();
       ++trackIter)
  {
    /// Get the track information
    const unsigned int trackKey = trackIter->first;
    const Trajectory traj = trackIter->second;
    SvtxTrackMap::Iter trackIter = m_trackMap->find(trackKey);
    SvtxTrack *track = trackIter->second;
    PHG4Particle *g4particle = trackeval->max_truth_particle_by_nclusters(track);
    ActsTrack actsProtoTrack = m_actsProtoTrackMap->find(trackKey)->second;
    
    const auto &[trackTips, mj] = traj.trajectory();
    m_trajNr = iTrack;

    /// Skip failed tracks
    if (trackTips.empty())
    {
      if (Verbosity() > 1)
        std::cout << "TrackTips empty in ActsEvaluator" << std::endl;
      continue;
    }

    if (trackTips.size() > 1)
    {
      std::cout << "There should not be a track with fit track tip > 1... Bailing."
                << std::endl;
      break;
    }

    auto &trackTip = trackTips.front();
    auto trajState =
        Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

    m_nMeasurements = trajState.nMeasurements;
    m_nStates = trajState.nStates;

    fillG4Particle(g4particle);
    fillProtoTrack(actsProtoTrack, topNode);
    fillFittedTrackParams(traj);
    visitTrackStates(traj, topNode);

    m_trackTree->Fill();

    /// Start fresh for the next track
    clearTrackVariables();
    ++iTrack;
  }

  m_eventNr++;

  if (Verbosity() > 1)
    std::cout << "Finished ActsEvaluator::process_event" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int ActsEvaluator::End(PHCompositeNode *topNode)
{
  m_trackFile->cd();
  m_trackTree->Write();
  m_trackFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

int ActsEvaluator::ResetEvent(PHCompositeNode *topNode)
{
  m_trajNr = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}


void ActsEvaluator::visitTrackStates(const Trajectory traj, PHCompositeNode *topNode)
{

  const auto &[trackTips, mj] = traj.trajectory();
  auto &trackTip = trackTips.front();

  mj.visitBackwards(trackTip, [&](const auto &state) {
    /// Only fill the track states with non-outlier measurement
    auto typeFlags = state.typeFlags();
    if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag))
    {
      return true;
    }

    /// Get the geometry ID
    auto geoID = state.referenceSurface().geoID();
    m_volumeID.push_back(geoID.volume());
    m_layerID.push_back(geoID.layer());
    m_moduleID.push_back(geoID.sensitive());

    auto meas = std::get<Measurement>(*state.uncalibrated());

    /// Get local position
    Acts::Vector2D local(meas.parameters()[Acts::ParDef::eLOC_0],
                         meas.parameters()[Acts::ParDef::eLOC_1]);
    /// Get global position
    Acts::Vector3D global(0, 0, 0);
    /// This is an arbitrary vector. Doesn't matter in coordinate transformation
    /// in Acts code
    Acts::Vector3D mom(1, 1, 1);
    meas.referenceSurface().localToGlobal(m_tGeometry->geoContext,
                                          local, mom, global);

    /// Get measurement covariance
    auto cov = meas.covariance();

    m_lx_hit.push_back(local.x());
    m_ly_hit.push_back(local.y());
    m_x_hit.push_back(global.x());
    m_y_hit.push_back(global.y());
    m_z_hit.push_back(global.z());

    /// Get the truth hit corresponding to this trackState
    /// We go backwards from hitID -> TrkrDefs::cluskey to g4hit with
    /// the map created in PHActsSourceLinks
    const unsigned int hitId = state.uncalibrated().hitID();
    float gt = -9999;
    Acts::Vector3D globalTruthPos = getGlobalTruthHit(topNode, hitId, gt);
    float gx = globalTruthPos(0);
    float gy = globalTruthPos(1);
    float gz = globalTruthPos(2);

    /// Get local truth position
    Acts::Vector2D truthlocal;

    const float r = sqrt(gx * gx + gy * gy + gz * gz);
    Acts::Vector3D globalTruthUnitDir(gx / r, gy / r, gz / r);

    meas.referenceSurface().globalToLocal(
        m_tGeometry->geoContext,
        globalTruthPos,
        globalTruthUnitDir,
        truthlocal);

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

    truthLOC0 = truthlocal.x();
    truthLOC1 = truthlocal.y();
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
      Acts::BoundParameters parameter(
          m_tGeometry->geoContext,
          state.predictedCovariance(),
          state.predicted(),
          state.referenceSurface().getSharedPtr());
      auto covariance = state.predictedCovariance();

      /// Local hit residual info
      auto H = meas.projector();
      auto resCov = cov + H * covariance * H.transpose();
      auto residual = meas.residual(parameter);
      m_res_x_hit.push_back(residual(Acts::ParDef::eLOC_0));
      m_res_y_hit.push_back(residual(Acts::ParDef::eLOC_1));
      m_err_x_hit.push_back(
          sqrt(resCov(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_err_y_hit.push_back(
          sqrt(resCov(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_pull_x_hit.push_back(
          residual(Acts::ParDef::eLOC_0) /
          sqrt(resCov(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_pull_y_hit.push_back(
          residual(Acts::ParDef::eLOC_1) /
          sqrt(resCov(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_dim_hit.push_back(state.calibratedSize());

      /// Predicted parameter
      m_eLOC0_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]);
      m_eLOC1_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
      m_ePHI_prt.push_back(parameter.parameters()[Acts::ParDef::ePHI]);
      m_eTHETA_prt.push_back(parameter.parameters()[Acts::ParDef::eTHETA]);
      m_eQOP_prt.push_back(parameter.parameters()[Acts::ParDef::eQOP]);
      m_eT_prt.push_back(parameter.parameters()[Acts::ParDef::eT]);

      /// Predicted residual
      m_res_eLOC0_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] -
                                truthLOC0);
      m_res_eLOC1_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] -
                                truthLOC1);
      m_res_ePHI_prt.push_back(parameter.parameters()[Acts::ParDef::ePHI] -
                               truthPHI);
      m_res_eTHETA_prt.push_back(
          parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA);
      m_res_eQOP_prt.push_back(parameter.parameters()[Acts::ParDef::eQOP] -
                               truthQOP);
      m_res_eT_prt.push_back(parameter.parameters()[Acts::ParDef::eT] -
                             truthTIME);

      /// Predicted parameter Uncertainties
      m_err_eLOC0_prt.push_back(
          sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_err_eLOC1_prt.push_back(
          sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_err_ePHI_prt.push_back(
          sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
      m_err_eTHETA_prt.push_back(
          sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
      m_err_eQOP_prt.push_back(
          sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
      m_err_eT_prt.push_back(
          sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

      /// Predicted parameter pulls
      m_pull_eLOC0_prt.push_back(
          (parameter.parameters()[Acts::ParDef::eLOC_0] - truthLOC0) /
          sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_pull_eLOC1_prt.push_back(
          (parameter.parameters()[Acts::ParDef::eLOC_1] - truthLOC1) /
          sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_pull_ePHI_prt.push_back(
          (parameter.parameters()[Acts::ParDef::ePHI] - truthPHI) /
          sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
      m_pull_eTHETA_prt.push_back(
          (parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA) /
          sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
      m_pull_eQOP_prt.push_back(
          (parameter.parameters()[Acts::ParDef::eQOP] - truthQOP) /
          sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
      m_pull_eT_prt.push_back(
          (parameter.parameters()[Acts::ParDef::eT] - truthTIME) /
          sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

      m_x_prt.push_back(parameter.position().x());
      m_y_prt.push_back(parameter.position().y());
      m_z_prt.push_back(parameter.position().z());
      m_px_prt.push_back(parameter.momentum().x());
      m_py_prt.push_back(parameter.momentum().y());
      m_pz_prt.push_back(parameter.momentum().z());
      m_pT_prt.push_back(parameter.pT());
      m_eta_prt.push_back(eta(parameter.position()));
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
      Acts::BoundParameters parameter(
          m_tGeometry->geoContext,
          state.filteredCovariance(), state.filtered(),
          state.referenceSurface().getSharedPtr());
      auto covariance = state.filteredCovariance();

      m_eLOC0_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]);
      m_eLOC1_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
      m_ePHI_flt.push_back(parameter.parameters()[Acts::ParDef::ePHI]);
      m_eTHETA_flt.push_back(parameter.parameters()[Acts::ParDef::eTHETA]);
      m_eQOP_flt.push_back(parameter.parameters()[Acts::ParDef::eQOP]);
      m_eT_flt.push_back(parameter.parameters()[Acts::ParDef::eT]);

      m_res_eLOC0_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] -
                                truthLOC0);
      m_res_eLOC1_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] -
                                truthLOC1);
      m_res_ePHI_flt.push_back(parameter.parameters()[Acts::ParDef::ePHI] -
                               truthPHI);
      m_res_eTHETA_flt.push_back(parameter.parameters()[Acts::ParDef::eTHETA] -
                                 truthTHETA);
      m_res_eQOP_flt.push_back(parameter.parameters()[Acts::ParDef::eQOP] -
                               truthQOP);
      m_res_eT_flt.push_back(parameter.parameters()[Acts::ParDef::eT] -
                             truthTIME);

      /// Filtered parameter uncertainties
      m_err_eLOC0_flt.push_back(
          sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_err_eLOC1_flt.push_back(
          sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_err_ePHI_flt.push_back(
          sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
      m_err_eTHETA_flt.push_back(
          sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
      m_err_eQOP_flt.push_back(
          sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
      m_err_eT_flt.push_back(
          sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

      /// Filtered parameter pulls
      m_pull_eLOC0_flt.push_back(
          (parameter.parameters()[Acts::ParDef::eLOC_0] - truthLOC0) /
          sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_pull_eLOC1_flt.push_back(
          (parameter.parameters()[Acts::ParDef::eLOC_1] - truthLOC1) /
          sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_pull_ePHI_flt.push_back(
          (parameter.parameters()[Acts::ParDef::ePHI] - truthPHI) /
          sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
      m_pull_eTHETA_flt.push_back(
          (parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA) /
          sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
      m_pull_eQOP_flt.push_back(
          (parameter.parameters()[Acts::ParDef::eQOP] - truthQOP) /
          sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
      m_pull_eT_flt.push_back(
          (parameter.parameters()[Acts::ParDef::eT] - truthTIME) /
          sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

      /// Other filtered parameter info
      m_x_flt.push_back(parameter.position().x());
      m_y_flt.push_back(parameter.position().y());
      m_z_flt.push_back(parameter.position().z());
      m_px_flt.push_back(parameter.momentum().x());
      m_py_flt.push_back(parameter.momentum().y());
      m_pz_flt.push_back(parameter.momentum().z());
      m_pT_flt.push_back(parameter.pT());
      m_eta_flt.push_back(eta(parameter.position()));
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
      Acts::BoundParameters parameter(
          m_tGeometry->geoContext,
          state.smoothedCovariance(), state.smoothed(),
          state.referenceSurface().getSharedPtr());
      auto covariance = state.smoothedCovariance();

      m_eLOC0_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]);
      m_eLOC1_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
      m_ePHI_smt.push_back(parameter.parameters()[Acts::ParDef::ePHI]);
      m_eTHETA_smt.push_back(parameter.parameters()[Acts::ParDef::eTHETA]);
      m_eQOP_smt.push_back(parameter.parameters()[Acts::ParDef::eQOP]);
      m_eT_smt.push_back(parameter.parameters()[Acts::ParDef::eT]);

      m_res_eLOC0_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] -
                                truthLOC0);
      m_res_eLOC1_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] -
                                truthLOC1);
      m_res_ePHI_smt.push_back(parameter.parameters()[Acts::ParDef::ePHI] -
                               truthPHI);
      m_res_eTHETA_smt.push_back(parameter.parameters()[Acts::ParDef::eTHETA] -
                                 truthTHETA);
      m_res_eQOP_smt.push_back(parameter.parameters()[Acts::ParDef::eQOP] -
                               truthQOP);
      m_res_eT_smt.push_back(parameter.parameters()[Acts::ParDef::eT] -
                             truthTIME);

      /// Smoothed parameter uncertainties
      m_err_eLOC0_smt.push_back(
          sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_err_eLOC1_smt.push_back(
          sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_err_ePHI_smt.push_back(
          sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
      m_err_eTHETA_smt.push_back(
          sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
      m_err_eQOP_smt.push_back(
          sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
      m_err_eT_smt.push_back(
          sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

      /// Smoothed parameter pulls
      m_pull_eLOC0_smt.push_back(
          (parameter.parameters()[Acts::ParDef::eLOC_0] - truthLOC0) /
          sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0)));
      m_pull_eLOC1_smt.push_back(
          (parameter.parameters()[Acts::ParDef::eLOC_1] - truthLOC1) /
          sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1)));
      m_pull_ePHI_smt.push_back(
          (parameter.parameters()[Acts::ParDef::ePHI] - truthPHI) /
          sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI)));
      m_pull_eTHETA_smt.push_back(
          (parameter.parameters()[Acts::ParDef::eTHETA] - truthTHETA) /
          sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA)));
      m_pull_eQOP_smt.push_back(
          (parameter.parameters()[Acts::ParDef::eQOP] - truthQOP) /
          sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP)));
      m_pull_eT_smt.push_back(
          (parameter.parameters()[Acts::ParDef::eT] - truthTIME) /
          sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT)));

      m_x_smt.push_back(parameter.position().x());
      m_y_smt.push_back(parameter.position().y());
      m_z_smt.push_back(parameter.position().z());
      m_px_smt.push_back(parameter.momentum().x());
      m_py_smt.push_back(parameter.momentum().y());
      m_pz_smt.push_back(parameter.momentum().z());
      m_pT_smt.push_back(parameter.pT());
      m_eta_smt.push_back(eta(parameter.position()));
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

    return true;
  }   /// Finish lambda function
  );  /// Finish multi trajectory visitBackwards call

  return;
}
TrkrDefs::cluskey ActsEvaluator::getClusKey(const unsigned int hitID)
{
  TrkrDefs::cluskey clusKey = 0;
  /// Unfortunately the map is backwards for looking up cluster key from
  /// hit ID. So we need to iterate over it. There won't be duplicates since
  /// the cluster key and hit id are a one-to-one map
  std::map<TrkrDefs::cluskey, unsigned int>::iterator
      hitIter = m_hitIdClusKey->begin();
  while (hitIter != m_hitIdClusKey->end())
  {
    if (hitIter->second == hitID)
    {
      clusKey = hitIter->first;
      break;
    }
    ++hitIter;
  }

  return clusKey;
}


Acts::Vector3D ActsEvaluator::getGlobalTruthHit(PHCompositeNode *topNode, 
						const unsigned int hitID,
						float &_gt)
{
  SvtxClusterEval *clustereval = m_svtxEvalStack->get_cluster_eval();

  TrkrDefs::cluskey clusKey = getClusKey(hitID);
  
  const PHG4Hit *g4hit = clustereval->max_truth_hit_by_energy(clusKey);
  
  float layer = (float) TrkrDefs::getLayer(clusKey);
  float gx = -9999;
  float gy = -9999;
  float gz = -9999;
  float gt = -9999;
  float gedep = -9999;
  
  if (g4hit)
    {
      /// Cluster the associated truth hits within the same layer to get
      /// the truth cluster position
      std::set<PHG4Hit *> truth_hits = clustereval->all_truth_hits(clusKey);
      std::vector<PHG4Hit *> contributing_hits;
      std::vector<double> contributing_hits_energy;
      std::vector<std::vector<double>> contributing_hits_entry;
      std::vector<std::vector<double>> contributing_hits_exit;
      m_svtxEvaluator->LayerClusterG4Hits(topNode, truth_hits, 
					  contributing_hits,
                                          contributing_hits_energy, 
					  contributing_hits_entry,
                                          contributing_hits_exit, 
					  layer, gx, gy, gz, gt,
                                          gedep);
    }
  
  /// Convert to acts units of mm
  gx *= Acts::UnitConstants::cm;
  gy *= Acts::UnitConstants::cm;
  gz *= Acts::UnitConstants::cm;
  
  Acts::Vector3D globalPos(gx, gy, gz);
  _gt = gt;
  return globalPos;
  
}

void ActsEvaluator::fillProtoTrack(ActsTrack track, PHCompositeNode *topNode)
{
  FW::TrackParameters params = track.getTrackParams();
  std::vector<SourceLink> sourceLinks = track.getSourceLinks();
  
  Acts::Vector3D position = params.position();
  Acts::Vector3D momentum = params.momentum();
  m_protoTrackPx = momentum(0);
  m_protoTrackPy = momentum(1);
  m_protoTrackPz = momentum(2);
  m_protoTrackX  = position(0);
  m_protoTrackY  = position(1);
  m_protoTrackZ  = position(2);

  for(int i = 0; i < sourceLinks.size(); ++i)
    {
      /// Get source link global position
      Acts::Vector2D loc(sourceLinks.at(i).location()(0),
			 sourceLinks.at(i).location()(1));
      Acts::Vector3D globalPos(0,0,0);
      Acts::Vector3D mom(0,0,0);
      
      sourceLinks.at(i).referenceSurface().localToGlobal(
                                            m_tGeometry->geoContext,
					    loc,
					    mom, globalPos);

      m_SLx.push_back(globalPos(0));
      m_SLy.push_back(globalPos(1));
      m_SLz.push_back(globalPos(2));
      m_SL_lx.push_back(loc(0));
      m_SL_ly.push_back(loc(1));
      
      /// Get corresponding truth hit position
      const unsigned int hitID = sourceLinks.at(i).hitID();
      float gt = -9999;
      Acts::Vector3D globalTruthPos = getGlobalTruthHit(topNode, hitID, gt);
      float gx = globalTruthPos(0);
      float gy = globalTruthPos(1);
      float gz = globalTruthPos(2);
      
      /// Get local truth position
      Acts::Vector2D truthlocal;
      
      const float r = sqrt(gx * gx + gy * gy + gz * gz);
      Acts::Vector3D globalTruthUnitDir(gx / r, gy / r, gz / r);
      
      sourceLinks.at(i).referenceSurface().globalToLocal(
					    m_tGeometry->geoContext,
					    globalTruthPos,
					    globalTruthUnitDir,
					    truthlocal);
      m_t_SL_lx.push_back(truthlocal(0));
      m_t_SL_ly.push_back(truthlocal(1));
      m_t_SL_gx.push_back(gx);
      m_t_SL_gy.push_back(gy);
      m_t_SL_gz.push_back(gz);
      


    }

}

void ActsEvaluator::fillFittedTrackParams(const Trajectory traj)
{
  m_hasFittedParams = false;
  const auto &[trackTips, mj] = traj.trajectory();
  auto &trackTip = trackTips.front();

  /// If it has track parameters, fill the values
  if (traj.hasTrackParameters(trackTip))
  {
    m_hasFittedParams = true;
    const auto &boundParam = traj.trackParameters(trackTip);
    const auto &parameter = boundParam.parameters();
    const auto &covariance = *boundParam.covariance();
    m_eLOC0_fit = parameter[Acts::ParDef::eLOC_0];
    m_eLOC1_fit = parameter[Acts::ParDef::eLOC_1];
    m_ePHI_fit = parameter[Acts::ParDef::ePHI];
    m_eTHETA_fit = parameter[Acts::ParDef::eTHETA];
    m_eQOP_fit = parameter[Acts::ParDef::eQOP];
    m_eT_fit = parameter[Acts::ParDef::eT];
    m_err_eLOC0_fit =
        sqrt(covariance(Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_0));
    m_err_eLOC1_fit =
        sqrt(covariance(Acts::ParDef::eLOC_1, Acts::ParDef::eLOC_1));
    m_err_ePHI_fit = sqrt(covariance(Acts::ParDef::ePHI, Acts::ParDef::ePHI));
    m_err_eTHETA_fit =
        sqrt(covariance(Acts::ParDef::eTHETA, Acts::ParDef::eTHETA));
    m_err_eQOP_fit = sqrt(covariance(Acts::ParDef::eQOP, Acts::ParDef::eQOP));
    m_err_eT_fit = sqrt(covariance(Acts::ParDef::eT, Acts::ParDef::eT));

    m_px_fit = boundParam.momentum()(0);
    m_py_fit = boundParam.momentum()(1);
    m_pz_fit = boundParam.momentum()(2);
    m_x_fit  = boundParam.position()(0);
    m_y_fit  = boundParam.position()(1);
    m_z_fit  = boundParam.position()(2);

    return;
  }

  /// Otherwise mark it as a bad fit
  m_eLOC0_fit = -9999;
  m_eLOC1_fit = -9999;
  m_ePHI_fit = -9999;
  m_eTHETA_fit = -9999;
  m_eQOP_fit = -9999;
  m_eT_fit = -9999;
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

  return;
}
void ActsEvaluator::fillG4Particle(PHG4Particle *part)
{
  SvtxTruthEval *trutheval = m_svtxEvalStack->get_truth_eval();

  if (part)
  {
    m_t_barcode = part->get_track_id();
    const auto pid = part->get_pid();
    m_t_charge = pid < 0 ? -1 : 1;
    const auto vtx = trutheval->get_vertex(part);
    m_t_vx = vtx->get_x() * Acts::UnitConstants::cm;
    m_t_vy = vtx->get_y() * Acts::UnitConstants::cm;
    m_t_vz = vtx->get_z() * Acts::UnitConstants::cm;
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

int ActsEvaluator::getNodes(PHCompositeNode *topNode)
{
  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!m_truthInfo)
  {
    std::cout << PHWHERE << "PHG4TruthInfoContainer not found, cannot continue!"
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_hitIdClusKey = findNode::getClass<std::map<TrkrDefs::cluskey, unsigned int>>(topNode, "HitIDClusIDActsMap");

  if (!m_hitIdClusKey)
  {
    std::cout << PHWHERE << "No HitID:ClusKey map on node tree. Bailing."
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");

  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No Acts Tracking geometry on node tree. Bailing"
              << std::endl;

    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_actsFitResults = findNode::getClass<std::map<const unsigned int, Trajectory>>(topNode, "ActsFitResults");

  if (!m_actsFitResults)
  {
    std::cout << PHWHERE << "No Acts fit results on node tree. Bailing"
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

  m_actsProtoTrackMap = findNode::getClass<std::map<unsigned int, ActsTrack>>(topNode, "ActsTrackMap");
  if (!m_actsProtoTrackMap)
    {
      std::cout << PHWHERE << "No Acts proto tracks on node tree. Bailing."
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
  m_protoTrackX  = -9999.;
  m_protoTrackY  = -9999.;
  m_protoTrackZ  = -9999.;


  return;
}

void ActsEvaluator::initializeTree()
{
  m_trackFile = new TFile(Name().c_str(), "RECREATE");

  m_trackTree = new TTree("tracktree", "A tree with Acts KF track information");

  m_trackTree->Branch("event_nr", &m_eventNr);
  m_trackTree->Branch("traj_nr", &m_trajNr);
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


  m_trackTree->Branch("g_protoTrackPx", &m_protoTrackPx, "m_protoTrackPx/F");
  m_trackTree->Branch("g_protoTrackPy", &m_protoTrackPy, "m_protoTrackPy/F");
  m_trackTree->Branch("g_protoTrackPz", &m_protoTrackPz, "m_protoTrackPz/F");
  m_trackTree->Branch("g_protoTrackX", &m_protoTrackX, "m_protoTrackX/F");
  m_trackTree->Branch("g_protoTrackY", &m_protoTrackY, "m_protoTrackY/F");
  m_trackTree->Branch("g_protoTrackZ", &m_protoTrackZ, "m_protoTrackZ/F");
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

  m_trackTree->Branch("hasFittedParams", &m_hasFittedParams);
  m_trackTree->Branch("eLOC0_fit", &m_eLOC0_fit);
  m_trackTree->Branch("eLOC1_fit", &m_eLOC1_fit);
  m_trackTree->Branch("ePHI_fit", &m_ePHI_fit);
  m_trackTree->Branch("eTHETA_fit", &m_eTHETA_fit);
  m_trackTree->Branch("eQOP_fit", &m_eQOP_fit);
  m_trackTree->Branch("eT_fit", &m_eT_fit);
  m_trackTree->Branch("err_eLOC0_fit", &m_err_eLOC0_fit);
  m_trackTree->Branch("err_eLOC1_fit", &m_err_eLOC1_fit);
  m_trackTree->Branch("err_ePHI_fit", &m_err_ePHI_fit);
  m_trackTree->Branch("err_eTHETA_fit", &m_err_eTHETA_fit);
  m_trackTree->Branch("err_eQOP_fit", &m_err_eQOP_fit);
  m_trackTree->Branch("err_eT_fit", &m_err_eT_fit);
  m_trackTree->Branch("g_px_fit", &m_px_fit);
  m_trackTree->Branch("g_py_fit", &m_py_fit);
  m_trackTree->Branch("g_pz_fit", &m_pz_fit);
  m_trackTree->Branch("g_x_fit" , &m_x_fit);
  m_trackTree->Branch("g_y_fit" , &m_y_fit);
  m_trackTree->Branch("g_z_fit" , &m_z_fit);


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
