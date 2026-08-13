#include "TpcPolyTrackSeedConverter.h"

#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeedHelper.h>
#include <trackbase_historic/TrackSeed_v2.h>

#include <cmath>
#include <iostream>
#include <memory>

TpcPolyTrackSeedConverter::TpcPolyTrackSeedConverter(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("TPC_POLYTRACKS")
  , m_outputNodeName("TpcTrackSeedContainer")
{
}

int TpcPolyTrackSeedConverter::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) { return Fun4AllReturnCodes::ABORTRUN;
}
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) { return Fun4AllReturnCodes::ABORTRUN;
}
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyTrackSeedConverter::getNodes(PHCompositeNode* topNode)
{
  m_polyTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_inputNodeName);
  if (!m_polyTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions)
  {
    std::cerr << Name() << "::getNodes - missing " << m_crossingDecisionNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_geometry)
  {
    std::cerr << Name() << "::getNodes - missing ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyTrackSeedConverter::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_trackSeeds = findNode::getClass<TrackSeedContainer>(topNode, m_outputNodeName);
  if (!m_trackSeeds)
  {
    m_trackSeeds = new TrackSeedContainer_v1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_trackSeeds, m_outputNodeName, "PHObject");
    svtxNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool TpcPolyTrackSeedConverter::isValidTrack(const Tpc_PolyTrack* track) const
{
  if (!track || track->get_fit_status() == 0) { return false;
}
  if (track->size_cluster_keys() == 0) { return false;
}
  if (track->get_nclusters() != track->size_cluster_keys()) { return false;
}
  if (track->get_nclusters() <= 20) { return false;
}

  const double pt = std::hypot(track->get_px(), track->get_py());
  if (!std::isfinite(pt) || pt <= 0.1) { return false;
}

  for (unsigned int i = 0; i < track->size_cluster_keys(); ++i)
  {
    if (track->get_cluster_key(i) == TrkrDefs::CLUSKEYMAX) { return false;
}
  }

  return std::isfinite(track->get_seed_x0()) &&
         std::isfinite(track->get_seed_y0()) &&
         std::isfinite(track->get_seed_z0()) &&
         std::isfinite(track->get_seed_phi()) &&
         std::isfinite(track->get_seed_slope()) &&
         std::isfinite(track->get_seed_q_over_r()) &&
         std::fabs(track->get_seed_q_over_r()) > 0.0;
}

bool TpcPolyTrackSeedConverter::publishSeed(const Tpc_PolyTrack* track) const
{
  if (!track || !m_crossingDecisions || !m_geometry || !m_trackSeeds) { return false;
}

  const TpcCrossingDecision* crossing_decision =
      m_crossingDecisions->get_decision(track->get_source_assembled_track_id());
  if (!crossing_decision) { return false;
}

  const short crossing = crossing_decision->get_selected_crossing();
  const double drift_velocity = m_geometry->get_drift_velocity();
  if (drift_velocity <= 0.0 || !std::isfinite(drift_velocity)) { return false;
}

  unsigned int side = 2;
  for (unsigned int i = 0; i < track->size_cluster_keys(); ++i)
  {
    const TrkrDefs::cluskey cluster_key = track->get_cluster_key(i);
    if (TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::tpcId) { continue;
}
    side = TpcDefs::getSide(cluster_key);
    break;
  }
  if (side > 1) { return false;
}

  const double z_bunch_separation = m_crossingPeriodNs * drift_velocity;
  const double z0_uncorrected = (side == 0) ?
      track->get_seed_z0() + static_cast<double>(crossing) * z_bunch_separation :
      track->get_seed_z0() - static_cast<double>(crossing) * z_bunch_separation;
  if (!std::isfinite(z0_uncorrected)) { return false;
}

  const double q_over_r = track->get_seed_q_over_r();
  const double helix_radius = std::fabs(1.0 / q_over_r);
  const double helix_sign = q_over_r > 0.0 ? 1.0 : -1.0;
  const double helix_center_x = track->get_seed_x0() - helix_sign * helix_radius * std::sin(track->get_seed_phi());
  const double helix_center_y = track->get_seed_y0() + helix_sign * helix_radius * std::cos(track->get_seed_phi());
  if (!std::isfinite(helix_center_x) || !std::isfinite(helix_center_y)) { return false;
}

  auto seed = std::make_unique<TrackSeed_v2>();
  seed->set_X0(static_cast<float>(helix_center_x));
  seed->set_Y0(static_cast<float>(helix_center_y));
  seed->set_Z0(static_cast<float>(z0_uncorrected));
  seed->set_phi(static_cast<float>(track->get_seed_phi()));
  seed->set_slope(static_cast<float>(track->get_seed_slope()));
  seed->set_qOverR(static_cast<float>(q_over_r));
  seed->set_crossing(crossing);

  for (unsigned int i = 0; i < track->size_cluster_keys(); ++i)
  {
    seed->insert_cluster_key(track->get_cluster_key(i));
  }

  const auto* converted_seed = m_trackSeeds->insert(seed.get());

  if (Verbosity() > 1 && converted_seed)
  {
    std::cout << Name() << "::publishSeed"
              << " track_id=" << track->get_track_id()
              << " source_track=" << track->get_source_assembled_track_id()
              << " crossing=" << crossing
              << " side=" << side
              << " drift_velocity=" << drift_velocity
              << " helix_radius=" << helix_radius
              << " poly_seed=(x0=" << track->get_seed_x0()
              << ", y0=" << track->get_seed_y0()
              << ", z0=" << track->get_seed_z0()
              << ", phi=" << track->get_seed_phi()
              << ", slope=" << track->get_seed_slope()
              << ", qOverR=" << track->get_seed_q_over_r() << ")"
              << " converted_seed=(x0=" << converted_seed->get_X0()
              << ", y0=" << converted_seed->get_Y0()
              << ", z0=" << converted_seed->get_Z0()
              << ", phi=" << converted_seed->get_phi()
              << ", slope=" << converted_seed->get_slope()
              << ", qOverR=" << converted_seed->get_qOverR() << ")"
              << " ncluster_keys=" << track->size_cluster_keys()
              << std::endl;
  }

  return converted_seed != nullptr;
}

int TpcPolyTrackSeedConverter::process_event(PHCompositeNode* topNode)
{
  if (!m_polyTracks || !m_trackSeeds || !m_crossingDecisions || !m_geometry)
  {
    if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  m_trackSeeds->Reset();

  unsigned int naccepted = 0;
  for (unsigned int itrack = 0; itrack < m_polyTracks->size(); ++itrack)
  {
    const Tpc_PolyTrack* track = m_polyTracks->get_track(itrack);
    if (!isValidTrack(track)) { continue;
}
    if (publishSeed(track)) { ++naccepted;
}
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - poly_tracks=" << m_polyTracks->size()
              << " tpc_seeds=" << naccepted << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
