#include "PHTpcResiduals.h"
#include "TpcSpaceChargeMatrixContainerv2.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrackState_v3.h>
#include <trackbase_historic/TrackSeed.h>

#include <trackreco/ActsPropagator.h>

#include <micromegas/MicromegasDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/GenericBoundTrackParameters.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Result.hpp>

#include <TFile.h>
#include <TH1.h>

#include <cassert>
#include <iostream>
#include <limits>

namespace
{

  // square
  template <class T>
  constexpr T square(const T& x)
  {
    return x * x;
  }

  // radius
  template <class T>
  T get_r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  template <class T>
  constexpr T deltaPhi(const T& phi)
  {
    if (phi > M_PI)
    {
      return phi - 2. * M_PI;
    }
    if (phi <= -M_PI)
    {
      return phi + 2. * M_PI;
    }

    return phi;
  }

  /// get sector median angle associated to a given index
  /** this assumes that sector 0 is centered on phi=0, then numbered along increasing phi */
  constexpr double get_sector_phi(int isec)
  {
    return isec * M_PI / 6;
  }

  // specify bins for which one will save histograms
  const std::vector<float> phi_rec = {get_sector_phi(9)};
  const std::vector<float> z_rec = {5.};

  //! get cluster keys from a given track
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }
    return out;
  }

  std::vector<TrkrDefs::cluskey> get_state_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (auto state_iter = track->begin_states();
         state_iter != track->end_states();
         ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      auto stateckey = tstate->get_cluskey();
      out.push_back(stateckey);
    }
    return out;
  }

  /// return number of clusters of a given type that belong to a tracks
  template <int type>
  int count_clusters(const std::vector<TrkrDefs::cluskey>& keys)
  {
    return std::count_if(keys.begin(), keys.end(),
                         [](const TrkrDefs::cluskey& key)
                         { return TrkrDefs::getTrkrId(key) == type; });
  }

  /// streamer
  std::ostream& operator<<(std::ostream& out, const Acts::Vector3& v)
  {
    out << "(" << v.x() << "," << v.y() << "," << v.z() << ")";
    return out;
  }

}  // namespace

//___________________________________________________________________________________
PHTpcResiduals::PHTpcResiduals(const std::string& name)
  : SubsysReco(name)
  , m_matrix_container(new TpcSpaceChargeMatrixContainerv2)
{
}

//___________________________________________________________________________________
int PHTpcResiduals::Init(PHCompositeNode* /*topNode*/)
{
  // configuration printout
  std::cout << "PHTpcResiduals::Init - m_maxTAlpha: " << m_maxTAlpha << std::endl;
  std::cout << "PHTpcResiduals::Init - m_maxTBeta: " << m_maxTBeta << std::endl;
  std::cout << "PHTpcResiduals::Init - m_maxResidualDrphi: " << m_maxResidualDrphi << " cm" << std::endl;
  std::cout << "PHTpcResiduals::Init - m_maxResidualDz: " << m_maxResidualDz << " cm" << std::endl;
  std::cout << "PHTpcResiduals::Init - m_minRPhiErr: " << m_minRPhiErr << " cm" << std::endl;
  std::cout << "PHTpcResiduals::Init - m_minZErr: " << m_minZErr << " cm" << std::endl;
  std::cout << "PHTpcResiduals::Init - m_minPt: " << m_minPt << " GeV/c" << std::endl;
  std::cout << "PHTpcResiduals::Init - m_requireCrossing: " << m_requireCrossing << std::endl;

  // reset counters
  m_total_tracks = 0;
  m_accepted_tracks = 0;

  m_total_clusters = 0;
  m_accepted_clusters = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHTpcResiduals::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  m_zMax = m_tGeometry->get_max_driftlength() + m_tGeometry->get_CM_halfwidth();
  m_zMin = -m_zMax;
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHTpcResiduals::process_event(PHCompositeNode* topNode)
{
  const auto returnVal = processTracks(topNode);
  ++m_event;

  return returnVal;
}

//___________________________________________________________________________________
int PHTpcResiduals::End(PHCompositeNode* /*topNode*/)
{
  std::cout << "PHTpcResiduals::End - writing matrices to " << m_outputfile << std::endl;

  // save matrix container in output file
  if (m_matrix_container)
  {
    std::unique_ptr<TFile> outputfile(TFile::Open(m_outputfile.c_str(), "RECREATE"));
    outputfile->cd();
    m_matrix_container->Write("TpcSpaceChargeMatrixContainer");
  }

  // print counters
  std::cout
      << "PHTpcResiduals::End -"
      << " track statistics total: " << m_total_tracks
      << " accepted: " << m_accepted_tracks
      << " fraction: " << 100. * m_accepted_tracks / m_total_tracks << "%"
      << std::endl;

  std::cout
      << "PHTpcResiduals::End -"
      << " cluster statistics total: " << m_total_clusters
      << " accepted: " << m_accepted_clusters << " fraction: "
      << 100. * m_accepted_clusters / m_total_clusters << "%"
      << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHTpcResiduals::processTracks(PHCompositeNode* /*topNode*/)
{
  if (Verbosity())
  {
    std::cout << "PHTpcResiduals::processTracks - proto track size " << m_trackMap->size() << std::endl;
  }

  for (const auto& [trackKey, track] : *m_trackMap)
  {
    if (Verbosity() > 1)
    {
      std::cout << "PHTpcResiduals::processTracks - Processing track key " << trackKey << std::endl;
    }

    ++m_total_tracks;
    if (checkTrack(track))
    {
      ++m_accepted_tracks;
      processTrack(track);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
bool PHTpcResiduals::checkTrack(SvtxTrack* track) const
{
  if (Verbosity() > 2)
  {
    std::cout << "PHTpcResiduals::checkTrack - pt: " << track->get_pt() << std::endl;
  }

  if (m_requireCrossing && track->get_crossing() != 0)
  {
    return false;
  }

  if (track->get_pt() < m_minPt)
  {
    return false;
  }

  // ignore tracks with too few mvtx, intt and micromegas hits
  const auto cluster_keys(get_cluster_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(cluster_keys) < 3)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(cluster_keys) < 2)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::tpcId>(cluster_keys) < 35)
  {
    return false;
  }
  //  if (m_useMicromegas && count_clusters<TrkrDefs::micromegasId>(cluster_keys) < 2)
  //  {
  //    return false;
  //  }

  const auto state_keys(get_state_keys(track));
  if (count_clusters<TrkrDefs::mvtxId>(state_keys) < 3)
  {
    return false;
  }
  if (count_clusters<TrkrDefs::inttId>(state_keys) < 2)
  {
    return false;
  }
  //  if (m_useMicromegas && count_clusters<TrkrDefs::micromegasId>(state_keys) < 2)
  //  {
  //    return false;
  //  }

  if (m_useMicromegas && checkTPOTResidual(track) == false)
  {
    return false;
  }

  return true;
}

//___________________________________________________________________________________
bool PHTpcResiduals::checkTPOTResidual(SvtxTrack* track) const
{
  bool flag = true;

  int nTPOTcluster = 0;
  int nTPOTstate = 0;
  int TPOTtileID = -1;
  for (const auto& cluskey : get_cluster_keys(track))
  {
    // make sure cluster is from TPOT
    const auto detId = TrkrDefs::getTrkrId(cluskey);
    if (detId != TrkrDefs::micromegasId)
    {
      continue;
    }
    TPOTtileID = MicromegasDefs::getTileId(cluskey);
    nTPOTcluster++;

    auto* const cluster = m_clusterContainer->findCluster(cluskey);

    SvtxTrackState* state = nullptr;

    // the track states from the Acts fit are fitted to fully corrected clusters, and are on the surface
    for (auto state_iter = track->begin_states();
         state_iter != track->end_states();
         ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      auto stateckey = tstate->get_cluskey();
      if (stateckey == cluskey)
      {
        state = tstate;
        break;
      }
    }

    const auto layer = TrkrDefs::getLayer(cluskey);

    if (!state)
    {
      if (Verbosity() > 1)
      {
        std::cout << "   no state for cluster " << cluskey << "  in layer " << layer << std::endl;
      }
      continue;
    }
    nTPOTstate++;

    const auto crossing = track->get_crossing();
    assert(crossing != std::numeric_limits<short>::max());

    // calculate residuals with respect to cluster
    // Get all the relevant information for residual calculation
    const auto globClusPos = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluskey, cluster, crossing);
    const double clusR = get_r(globClusPos(0), globClusPos(1));
    const double clusPhi = std::atan2(globClusPos(1), globClusPos(0));
    const double clusZ = globClusPos(2);

    const double globStateX = state->get_x();
    const double globStateY = state->get_y();
    const double globStateZ = state->get_z();
    const double globStatePx = state->get_px();
    const double globStatePy = state->get_py();
    const double globStatePz = state->get_pz();

    const double trackR = std::sqrt(square(globStateX) + square(globStateY));

    const double dr = clusR - trackR;
    const double trackDrDt = (globStateX * globStatePx + globStateY * globStatePy) / trackR;
    const double trackDxDr = globStatePx / trackDrDt;
    const double trackDyDr = globStatePy / trackDrDt;
    const double trackDzDr = globStatePz / trackDrDt;

    const double trackX = globStateX + dr * trackDxDr;
    const double trackY = globStateY + dr * trackDyDr;
    const double trackZ = globStateZ + dr * trackDzDr;
    const double trackPhi = std::atan2(trackY, trackX);

    // Calculate residuals
    // need to be calculated in local coordinates in the future
    const double drphi = clusR * deltaPhi(clusPhi - trackPhi);
    if (std::isnan(drphi))
    {
      continue;
    }

    const double dz = clusZ - trackZ;
    if (std::isnan(dz))
    {
      continue;
    }

    if (Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::checkTPOTResidual -"
                << " drphi: " << drphi
                << " dz: " << dz
                << std::endl;
    }

    // check rphi residual for layer 55
    if (layer == 55 && std::fabs(drphi) > 0.1)
    {
      flag = false;
      break;
    }

    // check z residual for layer 56
    if (layer == 56 && std::fabs(dz) > 1)
    {
      flag = false;
      break;
    }
  }

  if (flag)
  {
    // SCOZ has a half dead tile
    // only require one TPOT cluster/state from SCOP
    if (TPOTtileID == 0)
    {
      if (nTPOTcluster < 1 || nTPOTstate < 1)
      {
        flag = false;
      }
    }
    else if (TPOTtileID > 0)
    {
      if (nTPOTcluster < 2 || nTPOTstate < 2)
      {
        flag = false;
      }
    }
    else if (TPOTtileID < 0)
    {
      flag = false;
    }
  }

  return flag;
}

//___________________________________________________________________________________
Acts::BoundTrackParameters PHTpcResiduals::makeTrackParams(SvtxTrack* track) const
{
  return makeTrackParams(track, track->get_state(0));
}

//___________________________________________________________________________________
Acts::BoundTrackParameters PHTpcResiduals::makeTrackParams(SvtxTrack* track, SvtxTrackState* state) const
{
  Acts::Vector3 momentum(track->get_px(),
                         track->get_py(),
                         track->get_pz());
  double trackQ = track->get_charge() * Acts::UnitConstants::e;
  double p = state->get_p();

  // covariance
  const auto cov = m_transformer.rotateSvtxTrackCovToActs(state);

  // position
  const Acts::Vector3 position(
      state->get_x() * Acts::UnitConstants::cm,
      state->get_y() * Acts::UnitConstants::cm,
      state->get_z() * Acts::UnitConstants::cm);

  if (Verbosity())
  {
    std::cout << "PHTpcResiduals::makeTrackParams -"
              << " position: " << position
              << " momentum: " << momentum
              << std::endl;

    // covariance
    std::cout << "PHTpcResiduals::makeTrackParams - cov: " << std::endl;
    for (int i = 0; i < 6; ++i)
    {
      std::cout << "  ";
      for (int j = 0; j < 6; ++j)
      {
        std::cout << state->get_error(i, j) << " ";
      }
      std::cout << std::endl;
    }
  }

  const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(position);
  const auto actsFourPos = Acts::Vector4(position(0), position(1), position(2), 10 * Acts::UnitConstants::ns);
  return Acts::BoundTrackParameters::create(perigee, m_tGeometry->geometry().getGeoContext(),
                                            actsFourPos, momentum,
                                            trackQ / p, cov,
                                            Acts::ParticleHypothesis::pion())
      .value();
}

//_____________________________________________________________________________________________
void PHTpcResiduals::processTrack(SvtxTrack* track)
{
  if (Verbosity() > 1)
  {
    std::cout << "PHTpcResiduals::processTrack -"
              << " track momentum: " << track->get_p()
              << " position: " << Acts::Vector3(track->get_x(), track->get_y(), track->get_z())
              << std::endl;
  }
  ActsPropagator propagator(m_tGeometry);

  // create ACTS parameters from track parameters at origin
  auto trackParams = makeTrackParams(track);

  // store crossing. It is used in calculating cluster's global position
  const auto crossing = track->get_crossing();
  assert(crossing != std::numeric_limits<short>::max());

  for (const auto& cluskey : get_cluster_keys(track))
  {
    // increment counter
    ++m_total_clusters;

    // make sure cluster is from TPC
    const auto detId = TrkrDefs::getTrkrId(cluskey);
    if (detId != TrkrDefs::tpcId)
    {
      continue;
    }

    const auto layer = TrkrDefs::getLayer(cluskey);
    auto* const cluster = m_clusterContainer->findCluster(cluskey);
    const auto surface = m_tGeometry->maps().getSurface(cluskey, cluster);
    // auto result = propagator.propagateTrack(trackParams, surface);//surface aborter
    auto result = propagator.propagateTrack(trackParams, layer);  // layer aborter

    // skip if propagation failed
    if (!result.ok())
    {
      if (Verbosity() > 1)
      {
        std::cout << "Starting track params position/momentum: "
                  << trackParams.position(m_tGeometry->geometry().geoContext).transpose()
                  << std::endl
                  << std::endl
                  << "Track params phi/eta "
                  << std::atan2(trackParams.momentum().y(),
                                trackParams.momentum().x())
                  << " and "
                  << std::atanh(trackParams.momentum().z() / trackParams.momentum().norm())
                  << std::endl;
      }

      continue;
    }

    // get extrapolated track state, convert to sPHENIX and add to track
    auto& [pathLength, trackStateParams] = result.value();
    pathLength /= Acts::UnitConstants::cm;

    if (Verbosity() > 1)
    {
      std::cout << "PHTpcResiduals::processTrack -"
                << " path length: " << pathLength
                << " track momentum : "
                << trackParams.momentum()
                << " propagator momentum : "
                << trackStateParams.momentum()
                << std::endl;
    }

    addTrackState(track, cluskey, pathLength, trackStateParams);

    // calculate residuals with respect to cluster
    // Get all the relevant information for residual calculation
    const auto globClusPos = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluskey, cluster, crossing);
    const double clusR = get_r(globClusPos(0), globClusPos(1));
    const double clusPhi = std::atan2(globClusPos(1), globClusPos(0));
    const double clusZ = globClusPos(2);

    // cluster errors
    double clusRPhiErr = 0;
    double clusZErr = 0;

    clusRPhiErr = cluster->getRPhiError();
    clusZErr = cluster->getZError();

    if (Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack -"
                << " cluskey: " << cluskey
                << " clusR: " << clusR
                << " clusPhi: " << clusPhi << "+/-" << clusRPhiErr
                << " clusZ: " << clusZ << "+/-" << clusZErr
                << std::endl;
    }

    /*
     * as instructed by Christof, it should not be necessary to cut on small
     * cluster errors any more with clusters of version 4 or higher
     */

    const auto globalStatePos = trackStateParams.position(m_tGeometry->geometry().getGeoContext());
    const auto globalStateMom = trackStateParams.momentum();
    const auto globalStateCov = *trackStateParams.covariance();  // NOLINT(bugprone-unchecked-optional-access)

    const double trackRPhiErr = std::sqrt(globalStateCov(Acts::eBoundLoc0, Acts::eBoundLoc0)) / Acts::UnitConstants::cm;
    const double trackZErr = sqrt(globalStateCov(Acts::eBoundLoc1, Acts::eBoundLoc1)) / Acts::UnitConstants::cm;

    const double globStateX = globalStatePos.x() / Acts::UnitConstants::cm;
    const double globStateY = globalStatePos.y() / Acts::UnitConstants::cm;
    const double globStateZ = globalStatePos.z() / Acts::UnitConstants::cm;

    const double trackR = std::sqrt(square(globStateX) + square(globStateY));

    const double dr = clusR - trackR;
    const double trackDrDt = (globStateX * globalStateMom(0) + globStateY * globalStateMom(1)) / trackR;
    const double trackDxDr = globalStateMom(0) / trackDrDt;
    const double trackDyDr = globalStateMom(1) / trackDrDt;
    const double trackDzDr = globalStateMom(2) / trackDrDt;

    const double trackX = globStateX + dr * trackDxDr;
    const double trackY = globStateY + dr * trackDyDr;
    const double trackZ = globStateZ + dr * trackDzDr;
    const double trackPhi = std::atan2(trackY, trackX);

    if (Verbosity() > 2)
    {
      std::cout << "PHTpcResiduals::processTrack -"
                << " trackR: " << trackR
                << " dr: " << dr
                << " trackDrDt: " << trackDrDt
                << " trackDxDr: " << trackDxDr
                << " trackDyDr: " << trackDyDr
                << " trackDzDr: " << trackDzDr
                << " trackPhi: " << trackPhi << "+/-" << trackRPhiErr
                << " track position: (" << trackX << ", " << trackY << ", " << trackZ << ")"
                << std::endl;
    }

    const double erp = square(clusRPhiErr) + square(trackRPhiErr);
    if (std::isnan(erp))
    {
      continue;
    }

    const double ez = square(clusZErr) + square(trackZErr);
    if (std::isnan(ez))
    {
      continue;
    }

    // Calculate residuals
    const double drphi = clusR * deltaPhi(clusPhi - trackPhi);
    if (std::isnan(drphi))
    {
      continue;
    }

    const double dz = clusZ - trackZ;
    if (std::isnan(dz))
    {
      continue;
    }

    if (Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack -"
                << " drphi: " << drphi
                << " dz: " << dz
                << " erp: " << erp
                << " ez: " << ez
                << std::endl;
    }

    // check rphi and z error
    if (std::sqrt(erp) < m_minRPhiErr)
    {
      continue;
    }

    if (std::sqrt(ez) < m_minZErr)
    {
      continue;
    }

    const double trackPPhi = -trackStateParams.momentum()(0) * std::sin(trackPhi) + trackStateParams.momentum()(1) * std::cos(trackPhi);
    const double trackPR = trackStateParams.momentum()(0) * std::cos(trackPhi) + trackStateParams.momentum()(1) * std::sin(trackPhi);
    const double trackPZ = trackStateParams.momentum()(2);

    const double trackAlpha = -trackPPhi / trackPR;
    if (std::isnan(trackAlpha))
    {
      continue;
    }

    const double trackBeta = -trackPZ / trackPR;
    if (std::isnan(trackBeta))
    {
      continue;
    }

    if (Verbosity() > 3)
    {
      std::cout
          << "PHTpcResiduals::processTrack -"
          << " trackPPhi: " << trackPPhi
          << " trackPR: " << trackPR
          << " trackPZ: " << trackPZ
          << " trackAlpha: " << trackAlpha
          << " trackBeta: " << trackBeta
          << std::endl;
    }

    // get cell index
    const auto index = getCell(globClusPos);

    if (Verbosity() > 3)
    {
      std::cout << "Bin index found is " << index << std::endl;
    }

    if (index < 0)
    {
      continue;
    }

    if (Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack - layer: " << (int) TrkrDefs::getLayer(cluskey) << std::endl;
      std::cout << "PHTpcResiduals::processTrack -"
                << " cluster: (" << clusR << ", " << clusR * clusPhi << ", " << clusZ << ")"
                << " (" << clusRPhiErr << ", " << clusZErr << ")"
                << std::endl;

      std::cout << "PHTpcResiduals::processTrack -"
                << " track: (" << trackR << ", " << clusR * trackPhi << ", " << trackZ << ")"
                << " (" << trackAlpha << ", " << trackBeta << ")"
                << " (" << trackRPhiErr << ", " << trackZErr << ")"
                << std::endl;
      std::cout << std::endl;
    }

    // check track angles and residuals agains cuts
    if (std::abs(trackAlpha) > m_maxTAlpha || std::abs(drphi) > m_maxResidualDrphi)
    {
      continue;
    }

    if (std::abs(trackBeta) > m_maxTBeta || std::abs(dz) > m_maxResidualDz)
    {
      continue;
    }

    // Fill distortion matrices
    m_matrix_container->add_to_lhs(index, 0, 0, square(clusR) / erp);
    m_matrix_container->add_to_lhs(index, 0, 1, 0);
    m_matrix_container->add_to_lhs(index, 0, 2, clusR * trackAlpha / erp);

    m_matrix_container->add_to_lhs(index, 1, 0, 0);
    m_matrix_container->add_to_lhs(index, 1, 1, 1. / ez);
    m_matrix_container->add_to_lhs(index, 1, 2, trackBeta / ez);

    m_matrix_container->add_to_lhs(index, 2, 0, clusR * trackAlpha / erp);
    m_matrix_container->add_to_lhs(index, 2, 1, trackBeta / ez);
    m_matrix_container->add_to_lhs(index, 2, 2, square(trackAlpha) / erp + square(trackBeta) / ez);

    m_matrix_container->add_to_rhs(index, 0, clusR * drphi / erp);
    m_matrix_container->add_to_rhs(index, 1, dz / ez);
    m_matrix_container->add_to_rhs(index, 2, trackAlpha * drphi / erp + trackBeta * dz / ez);

    // also update rphi reduced matrices
    m_matrix_container->add_to_lhs_rphi(index, 0, 0, square(clusR) / erp);
    m_matrix_container->add_to_lhs_rphi(index, 0, 1, clusR * trackAlpha / erp);
    m_matrix_container->add_to_lhs_rphi(index, 1, 0, clusR * trackAlpha / erp);
    m_matrix_container->add_to_lhs_rphi(index, 1, 1, square(trackAlpha) / erp);

    m_matrix_container->add_to_rhs_rphi(index, 0, clusR * drphi / erp);
    m_matrix_container->add_to_rhs_rphi(index, 1, trackAlpha * drphi / erp);

    // also update z reduced matrices
    m_matrix_container->add_to_lhs_z(index, 0, 0, 1. / ez);
    m_matrix_container->add_to_lhs_z(index, 0, 1, trackBeta / ez);
    m_matrix_container->add_to_lhs_z(index, 1, 0, trackBeta / ez);
    m_matrix_container->add_to_lhs_z(index, 1, 1, square(trackBeta) / ez);

    m_matrix_container->add_to_rhs_z(index, 0, dz / ez);
    m_matrix_container->add_to_rhs_z(index, 1, trackBeta * dz / ez);

    // update entries in cell
    m_matrix_container->add_to_entries(index);

    // increment number of accepted clusters
    ++m_accepted_clusters;
  }
}

//_______________________________________________________________________________________________________
void PHTpcResiduals::addTrackState(SvtxTrack* track, TrkrDefs::cluskey key, float pathlength, const Acts::BoundTrackParameters& params)
{
  /* this is essentially a copy of the code from trackbase_historic/ActsTransformations::fillSvtxTrackStates */

  // create track state
  SvtxTrackState_v3 state(pathlength);

  // save global position
  const auto global = params.position(m_tGeometry->geometry().getGeoContext());
  state.set_x(global.x() / Acts::UnitConstants::cm);
  state.set_y(global.y() / Acts::UnitConstants::cm);
  state.set_z(global.z() / Acts::UnitConstants::cm);

  // save momentum
  const auto momentum = params.momentum();
  state.set_px(momentum.x());
  state.set_py(momentum.y());
  state.set_pz(momentum.z());

  // covariance
  const auto globalCov = m_transformer.rotateActsCovToSvtxTrack(params);
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      state.set_error(i, j, globalCov(i, j));
    }
  }

  state.set_name(std::to_string(key));
  state.set_cluskey(key);
  track->insert_state(&state);
}

//_______________________________________________________________________________
int PHTpcResiduals::getCell(const Acts::Vector3& loc)
{
  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions(phibins, rbins, zbins);

  // phi
  float phi = std::atan2(loc(1), loc(0));
  while (phi < m_phiMin)
  {
    phi += 2. * M_PI;
  }
  while (phi >= m_phiMax)
  {
    phi -= 2. * M_PI;
  }
  const int iphi = phibins * (phi - m_phiMin) / (m_phiMax - m_phiMin);

  // r
  const auto r = get_r(loc(0), loc(1));
  if (r < m_rMin || r >= m_rMax)
  {
    return -1;
  }
  const int ir = rbins * (r - m_rMin) / (m_rMax - m_rMin);

  // z
  const auto z = loc(2);
  if (z < m_zMin || z >= m_zMax)
  {
    return -1;
  }
  const int iz = zbins * (z - m_zMin) / (m_zMax - m_zMin);

  // get index from matrix container
  return m_matrix_container->get_cell_index(iphi, ir, iz);
}

//_______________________________________________________________________________
int PHTpcResiduals::createNodes(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________________________________________________
int PHTpcResiduals::getNodes(PHCompositeNode* topNode)
{
  // clusters
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterContainer)
  {
    std::cout << PHWHERE << "No TRKR_CLUSTER node on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << "ActsTrackingGeometry not on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tracks
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxSiliconMMTrackMap");
  if (!m_trackMap)
  {
    std::cout << PHWHERE << "SvtxSiliconMMTrackMap not on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);
  if (m_disable_module_edge_corr)
  {
    m_globalPositionWrapper.set_enable_module_edge_corr(false);
  }
  if (m_disable_static_corr)
  {
    m_globalPositionWrapper.set_enable_static_corr(false);
  }
  if (m_disable_average_corr)
  {
    m_globalPositionWrapper.set_enable_average_corr(false);
  }
  if (m_disable_fluctuation_corr)
  {
    m_globalPositionWrapper.set_enable_fluctuation_corr(false);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
void PHTpcResiduals::setGridDimensions(const int phiBins, const int rBins, const int zBins)
{
  m_matrix_container->set_grid_dimensions(phiBins, rBins, zBins);
}
