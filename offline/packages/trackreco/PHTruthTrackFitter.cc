#include "PHTruthTrackFitter.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrackState_v3.h>
#include <trackbase_historic/SvtxTrack_v4.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <utility>

namespace
{
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  bool is_finite(float value)
  {
    return std::isfinite(value);
  }

  float average_or(float a, float b, float fallback)
  {
    const bool aok = is_finite(a);
    const bool bok = is_finite(b);

    if (aok && bok)
    {
      return 0.5 * (a + b);
    }
    if (aok)
    {
      return a;
    }
    if (bok)
    {
      return b;
    }

    return fallback;
  }

  bool valid_track_id(unsigned int trackid)
  {
    return trackid != std::numeric_limits<unsigned int>::max();
  }
}  // namespace

PHTruthTrackFitter::PHTruthTrackFitter(const std::string& name)
  : SubsysReco(name)
{
}

int PHTruthTrackFitter::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "PHTruthTrackFitter::InitRun - output track map: " << m_trackMapName << std::endl;
  }

  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackFitter::process_event(PHCompositeNode* /*topNode*/)
{
  m_trackMap->Reset();

  unsigned int skipped_tracks = 0;
  for (auto seed : *m_seedMap)
  {
    if (!seed)
    {
      continue;
    }

    auto* tpc_seed = getSeed(m_tpcSeeds, seed->get_tpc_seed_index());
    auto* silicon_seed = getSeed(m_siliconSeeds, seed->get_silicon_seed_index());

    if (!tpc_seed && !silicon_seed)
    {
      ++skipped_tracks;
      continue;
    }

    const auto truth_track_id = getTruthTrackId(seed, tpc_seed, silicon_seed);
    if (!valid_track_id(truth_track_id))
    {
      if (Verbosity() > 1)
      {
        std::cout << "PHTruthTrackFitter::process_event - could not determine truth id for seed" << std::endl;
      }
      ++skipped_tracks;
      continue;
    }

    auto* g4particle = m_g4TruthInfo->GetParticle(truth_track_id);
    if (!g4particle)
    {
      if (Verbosity() > 1)
      {
        std::cout << "PHTruthTrackFitter::process_event - no PHG4Particle for track id "
                  << truth_track_id << std::endl;
      }
      ++skipped_tracks;
      continue;
    }

    const auto* g4vertex = m_g4TruthInfo->GetVtx(g4particle->get_vtx_id());
    if (!g4vertex)
    {
      if (Verbosity() > 1)
      {
        std::cout << "PHTruthTrackFitter::process_event - no PHG4VtxPoint for track id "
                  << truth_track_id << std::endl;
      }
      ++skipped_tracks;
      continue;
    }

    SvtxTrack_v4 track;
    track.set_tpc_seed(tpc_seed);
    track.set_silicon_seed(silicon_seed);
    track.set_crossing(getCrossing(tpc_seed, silicon_seed));
    track.set_vertex_id(g4particle->get_vtx_id());
    track.set_charge(getCharge(g4particle, tpc_seed, silicon_seed));
    track.set_chisq(0);

    track.set_x(g4vertex->get_x());
    track.set_y(g4vertex->get_y());
    track.set_z(g4vertex->get_z());
    track.set_px(g4particle->get_px());
    track.set_py(g4particle->get_py());
    track.set_pz(g4particle->get_pz());

    for (int i = 0; i < 6; ++i)
    {
      for (int j = i; j < 6; ++j)
      {
        track.set_error(i, j, 0);
      }
    }
    track.set_error(0, 0, square(m_positionError));
    track.set_error(1, 1, square(m_positionError));
    track.set_error(2, 2, square(m_zError));

    unsigned int state_index = 1;
    for (const auto* track_seed : {silicon_seed, tpc_seed})
    {
      if (!track_seed)
      {
        continue;
      }

      for (auto iter = track_seed->begin_cluster_keys(); iter != track_seed->end_cluster_keys(); ++iter)
      {
        if (addStateFromCluster(&track, *iter, truth_track_id, g4particle, g4vertex, state_index))
        {
          ++state_index;
        }
      }
    }

    if (track.size_states() <= 1)
    {
      if (Verbosity() > 1)
      {
        std::cout << "PHTruthTrackFitter::process_event - no truth states for track id "
                  << truth_track_id << std::endl;
      }
      ++skipped_tracks;
      continue;
    }

    track.set_ndf(std::max<int>(0, 2 * static_cast<int>(track.size_states()) - 5));

    const unsigned int track_id = m_trackMap->size();
    track.set_id(track_id);
    m_trackMap->insertWithKey(&track, track_id);
  }

  if (Verbosity() > 0)
  {
    std::cout << "PHTruthTrackFitter::process_event - built " << m_trackMap->size()
              << " truth tracks, skipped " << skipped_tracks << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackFitter::End(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackFitter::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  auto* dst_node = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    std::cerr << PHWHERE << "DST node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHNodeIterator dst_iter(dst_node);
  auto* svtx_node = dynamic_cast<PHCompositeNode*>(dst_iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtx_node)
  {
    svtx_node = new PHCompositeNode("SVTX");
    dst_node->addNode(svtx_node);
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    m_trackMap = new SvtxTrackMap_v2;
    auto* track_node = new PHIODataNode<PHObject>(m_trackMap, m_trackMapName, "PHObject");
    svtx_node->addNode(track_node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackFitter::getNodes(PHCompositeNode* topNode)
{
  m_seedMap = findNode::getClass<TrackSeedContainer>(topNode, m_svtxSeedMapName);
  if (!m_seedMap)
  {
    std::cout << PHWHERE << "No " << m_svtxSeedMapName << " on node tree. Bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!m_tpcSeeds)
  {
    std::cout << PHWHERE << "No TpcTrackSeedContainer on node tree. Bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!m_siliconSeeds)
  {
    std::cout << PHWHERE << "No SiliconTrackSeedContainer on node tree. Bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterMapName);
  if (!m_clusterMap)
  {
    std::cout << PHWHERE << "No " << m_clusterMapName << " on node tree. Bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_clusterHitMap = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterHitMap)
  {
    std::cout << PHWHERE << "No TRKR_CLUSTERHITASSOC on node tree. Bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_hitTruthAssoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!m_hitTruthAssoc)
  {
    std::cout << PHWHERE << "No TRKR_HITTRUTHASSOC on node tree. Bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_g4TruthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_g4TruthInfo)
  {
    std::cout << PHWHERE << "No G4TruthInfo on node tree. Bailing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_g4HitsTpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4HitsIntt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4HitsMvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4HitsMicromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  return Fun4AllReturnCodes::EVENT_OK;
}

TrackSeed* PHTruthTrackFitter::getSeed(TrackSeedContainer* container, unsigned int index) const
{
  if (!container || index >= container->size())
  {
    return nullptr;
  }

  return container->get(index);
}

unsigned int PHTruthTrackFitter::getTruthTrackId(const TrackSeed* svtxSeed,
                                                 const TrackSeed* tpcSeed,
                                                 const TrackSeed* siliconSeed) const
{
  auto truth_track_id = m_invalidTruthTrackId;
  for (const auto* seed : {svtxSeed, tpcSeed, siliconSeed})
  {
    if (!seed)
    {
      continue;
    }

    const auto seed_truth_track_id = seed->get_truth_track_id();
    if (!valid_track_id(seed_truth_track_id))
    {
      continue;
    }

    if (valid_track_id(truth_track_id) && seed_truth_track_id != truth_track_id)
    {
      if (Verbosity() > 0)
      {
        std::cout << "PHTruthTrackFitter::getTruthTrackId - inconsistent seed truth ids "
                  << truth_track_id << " and " << seed_truth_track_id << std::endl;
      }
      return m_invalidTruthTrackId;
    }

    truth_track_id = seed_truth_track_id;
  }

  return truth_track_id;
}

std::vector<const PHG4Hit*> PHTruthTrackFitter::getTruthHits(TrkrDefs::cluskey cluskey) const
{
  std::vector<const PHG4Hit*> truth_hits;
  if (!m_clusterHitMap || !m_hitTruthAssoc)
  {
    return truth_hits;
  }

  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);
  const auto trkrid = TrkrDefs::getTrkrId(hitsetkey);
  const auto hitrange = m_clusterHitMap->getHits(cluskey);

  std::set<PHG4HitDefs::keytype> used_g4hits;
  for (auto clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
  {
    const auto hitkey = clushititer->second;

    TrkrHitTruthAssoc::MMap temp_map;
    m_hitTruthAssoc->getG4Hits(hitsetkey, hitkey, temp_map);

    for (const auto& hit_truth_iter : temp_map)
    {
      const auto g4hitkey = hit_truth_iter.second.second;
      if (!used_g4hits.insert(g4hitkey).second)
      {
        continue;
      }

      const auto* g4hit = getG4Hit(trkrid, g4hitkey);
      if (g4hit)
      {
        truth_hits.push_back(g4hit);
      }
    }
  }

  return truth_hits;
}

const PHG4Hit* PHTruthTrackFitter::getG4Hit(unsigned int trkrid, PHG4HitDefs::keytype g4hitkey) const
{
  PHG4HitContainer* container = nullptr;
  switch (trkrid)
  {
  case TrkrDefs::tpcId:
    container = m_g4HitsTpc;
    break;
  case TrkrDefs::inttId:
    container = m_g4HitsIntt;
    break;
  case TrkrDefs::mvtxId:
    container = m_g4HitsMvtx;
    break;
  case TrkrDefs::micromegasId:
    container = m_g4HitsMicromegas;
    break;
  default:
    break;
  }

  return container ? container->findHit(g4hitkey) : nullptr;
}

bool PHTruthTrackFitter::addStateFromCluster(SvtxTrack* track,
                                             TrkrDefs::cluskey cluskey,
                                             unsigned int truthTrackId,
                                             const PHG4Particle* particle,
                                             const PHG4VtxPoint* vertex,
                                             unsigned int stateIndex) const
{
  if (!m_clusterMap->findCluster(cluskey))
  {
    return false;
  }

  double weight_sum = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double px = 0;
  double py = 0;
  double pz = 0;
  double local_x = 0;
  double local_y = 0;

  for (const auto* g4hit : getTruthHits(cluskey))
  {
    if (!g4hit || g4hit->get_trkid() != static_cast<int>(truthTrackId))
    {
      continue;
    }

    const auto hit_x = average_or(g4hit->get_x(0), g4hit->get_x(1), std::numeric_limits<float>::quiet_NaN());
    const auto hit_y = average_or(g4hit->get_y(0), g4hit->get_y(1), std::numeric_limits<float>::quiet_NaN());
    const auto hit_z = average_or(g4hit->get_z(0), g4hit->get_z(1), std::numeric_limits<float>::quiet_NaN());
    if (!is_finite(hit_x) || !is_finite(hit_y) || !is_finite(hit_z))
    {
      continue;
    }

    const auto hit_px = average_or(g4hit->get_px(0), g4hit->get_px(1), particle->get_px());
    const auto hit_py = average_or(g4hit->get_py(0), g4hit->get_py(1), particle->get_py());
    const auto hit_pz = average_or(g4hit->get_pz(0), g4hit->get_pz(1), particle->get_pz());
    const auto hit_local_x = average_or(g4hit->get_local_x(0), g4hit->get_local_x(1), 0);
    const auto hit_local_y = average_or(g4hit->get_local_y(0), g4hit->get_local_y(1), 0);

    double weight = g4hit->get_edep();
    if (!std::isfinite(weight) || weight <= 0)
    {
      weight = 1;
    }

    weight_sum += weight;
    x += weight * hit_x;
    y += weight * hit_y;
    z += weight * hit_z;
    px += weight * hit_px;
    py += weight * hit_py;
    pz += weight * hit_pz;
    local_x += weight * hit_local_x;
    local_y += weight * hit_local_y;
  }

  if (weight_sum <= 0)
  {
    return false;
  }

  x /= weight_sum;
  y /= weight_sum;
  z /= weight_sum;
  px /= weight_sum;
  py /= weight_sum;
  pz /= weight_sum;
  local_x /= weight_sum;
  local_y /= weight_sum;

  float pathlength = getPathLength(vertex, x, y, z, stateIndex);
  while (track->count_states(pathlength) != 0)
  {
    pathlength += 1.e-3;
  }

  SvtxTrackState_v3 state(pathlength);
  state.set_name("PHTruthTrackFitter");
  state.set_cluskey(cluskey);
  state.set_x(x);
  state.set_y(y);
  state.set_z(z);
  state.set_px(px);
  state.set_py(py);
  state.set_pz(pz);
  state.set_localX(local_x);
  state.set_localY(local_y);

  for (int i = 0; i < 6; ++i)
  {
    for (int j = i; j < 6; ++j)
    {
      state.set_error(i, j, 0);
    }
  }
  state.set_error(0, 0, square(m_positionError));
  state.set_error(1, 1, square(m_positionError));
  state.set_error(2, 2, square(m_zError));

  track->insert_state(&state);
  return true;
}

float PHTruthTrackFitter::getPathLength(const PHG4VtxPoint* vertex,
                                        float x, float y, float z,
                                        unsigned int stateIndex) const
{
  if (vertex)
  {
    const auto dx = x - vertex->get_x();
    const auto dy = y - vertex->get_y();
    const auto dz = z - vertex->get_z();
    const auto pathlength = std::sqrt(square(dx) + square(dy) + square(dz));
    if (std::isfinite(pathlength) && pathlength > 0)
    {
      return pathlength;
    }
  }

  return static_cast<float>(stateIndex);
}

int PHTruthTrackFitter::getCharge(const PHG4Particle* particle,
                                  const TrackSeed* tpcSeed,
                                  const TrackSeed* siliconSeed) const
{
  for (const auto* seed : {tpcSeed, siliconSeed})
  {
    if (!seed)
    {
      continue;
    }

    const auto charge = seed->get_charge();
    if (std::abs(charge) == 1)
    {
      return charge;
    }
  }

  const auto* pdg_particle = particle ? TDatabasePDG::Instance()->GetParticle(particle->get_pid()) : nullptr;
  if (pdg_particle && pdg_particle->Charge() < 0)
  {
    return -1;
  }

  return 1;
}

short int PHTruthTrackFitter::getCrossing(const TrackSeed* tpcSeed, const TrackSeed* siliconSeed) const
{
  for (const auto* seed : {siliconSeed, tpcSeed})
  {
    if (!seed)
    {
      continue;
    }

    const auto crossing = seed->get_crossing();
    if (crossing != std::numeric_limits<short int>::max())
    {
      return crossing;
    }
  }

  return m_defaultCrossing;
}
