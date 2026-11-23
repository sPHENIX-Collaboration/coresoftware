#include "KFParticle_detectorInfo.h"
#include "KFParticle_utilities.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>              // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>              // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>               // for getSectorId, getSide
#include <trackbase/TrkrCluster.h>           // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>  // for TrkrClusterContainer
#include <trackbase/TrkrDefs.h>              // for getLayer, getTrkrId

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, SvtxTrac...
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <phool/getClass.h>  // for getClass

#include <KFParticle.h>  // for KFParticle

#include <TTree.h>  // for TTree

#include <iostream>  // for operator<<, endl
#include <limits>
#include <map>
#include <string>

namespace
{
  std::map<std::string, int> Use =
      {
          {"MVTX", 1},
          {"INTT", 1},
          {"TPC", 1},
          {"TPOT", 1},
          {"EMCAL", 0},
          {"OHCAL", 0},
          {"IHCAL", 0}};
}  // namespace

void KFParticle_detectorInfo::initializeBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  if (m_config.detailedTracking)
  {
    m_tree->Branch((daughter_number + "_residual_x").c_str(), &residual_x[daughter_id]);
    m_tree->Branch((daughter_number + "_residual_y").c_str(), &residual_y[daughter_id]);
    m_tree->Branch((daughter_number + "_residual_z").c_str(), &residual_z[daughter_id]);
    m_tree->Branch((daughter_number + "_layer").c_str(), &detector_layer[daughter_id]);
  }

  for (const auto &subdetector : Use)
  {
    if (subdetector.second)
    {
      initializeSubDetectorBranches(m_tree, subdetector.first, daughter_id, daughter_number);
    }
  }
}
void KFParticle_detectorInfo::initializeSubDetectorBranches(TTree *m_tree, const std::string &detectorName, int daughter_id, const std::string &daughter_number)
{
  const std::string branchName = daughter_number + "_" + detectorName;
  if (detectorName == "MVTX")
  {
    if (m_config.detailedTracking)
    {
      m_tree->Branch((branchName + "_staveID").c_str(), &mvtx_staveID[daughter_id]);
      m_tree->Branch((branchName + "_chipID").c_str(), &mvtx_chipID[daughter_id]);
    }
    m_tree->Branch((branchName + "_nHits").c_str(), &detector_nHits_MVTX[daughter_id]);
    m_tree->Branch((branchName + "_nStates").c_str(), &detector_nStates_MVTX[daughter_id]);
  }
  else if (detectorName == "INTT")
  {
    if (m_config.detailedTracking)
    {
      m_tree->Branch((branchName + "_ladderZID").c_str(), &intt_ladderZID[daughter_id]);
      m_tree->Branch((branchName + "_ladderPhiID").c_str(), &intt_ladderPhiID[daughter_id]);
    }
    m_tree->Branch((branchName + "_nHits").c_str(), &detector_nHits_INTT[daughter_id]);
    m_tree->Branch((branchName + "_nStates").c_str(), &detector_nStates_INTT[daughter_id]);
  }
  else if (detectorName == "TPC")
  {
    if (m_config.detailedTracking)
    {
      m_tree->Branch((branchName + "_sectorID").c_str(), &tpc_sectorID[daughter_id]);
      m_tree->Branch((branchName + "_side").c_str(), &tpc_side[daughter_id]);
    }
    m_tree->Branch((branchName + "_nHits").c_str(), &detector_nHits_TPC[daughter_id]);
    m_tree->Branch((branchName + "_nStates").c_str(), &detector_nStates_TPC[daughter_id]);
  }
  else if (detectorName == "TPOT")
  {
    m_tree->Branch((branchName + "_nHits").c_str(), &detector_nHits_TPOT[daughter_id]);
    m_tree->Branch((branchName + "_nStates").c_str(), &detector_nStates_TPOT[daughter_id]);
  }
}

void KFParticle_detectorInfo::fill(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id)
{
  auto *dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_config.trackMapNodeName << " does not exist" << std::endl;
  }

  std::string geoName = "ActsGeometry";
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, geoName);
  if (!geometry)
  {
    std::cout << "KFParticle detector info: " << geoName << " does not exist" << std::endl;
  }

  auto *dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!dst_clustermap)
  {
    dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
    if (!dst_clustermap)
    {
      std::cout << "KFParticle detector info: TRKR_CLUSTER does not exist" << std::endl;
    }
  }

  SvtxTrack *track = KFParticle_utilities::getTrack(daughter.Id(), dst_trackmap);
  detector_nHits_MVTX[daughter_id] = 0;
  detector_nHits_INTT[daughter_id] = 0;
  detector_nHits_TPC[daughter_id] = 0;
  detector_nHits_TPOT[daughter_id] = 0;

  bool isTrack_v5 = track->InheritsFrom("SvtxTrack_v5");
  TrackSeedContainer *siliconSeedContainer = nullptr;
  TrackSeedContainer *tpcSeedContainer = nullptr;

  if (isTrack_v5)
  {
    siliconSeedContainer = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
    tpcSeedContainer = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  }

  TrackSeed *silseed = isTrack_v5 ? KFParticle_utilities::findSeed(siliconSeedContainer, track->get_silicon_seed_id()) : track->get_silicon_seed();
  TrackSeed *tpcseed = isTrack_v5 ? KFParticle_utilities::findSeed(tpcSeedContainer, track->get_tpc_seed_id()) : track->get_tpc_seed();

  if (silseed)
  {
    for (auto cluster_iter = silseed->begin_cluster_keys(); cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));

      unsigned int staveId;
      unsigned int chipId;
      unsigned int ladderZId;
      unsigned int ladderPhiId;
      unsigned int sectorId;
      unsigned int side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();

      if (Use["MVTX"] && trackerID == TrkrDefs::mvtxId)
      {
        staveId = MvtxDefs::getStaveId(cluster_key);
        chipId = MvtxDefs::getChipId(cluster_key);
        ++detector_nHits_MVTX[daughter_id];
      }
      else if (Use["INTT"] && trackerID == TrkrDefs::inttId)
      {
        ladderZId = InttDefs::getLadderZId(cluster_key);
        ladderPhiId = InttDefs::getLadderPhiId(cluster_key);
        ++detector_nHits_INTT[daughter_id];
      }

      if (m_config.detailedTracking)
      {
        mvtx_staveID[daughter_id].push_back(staveId);
        mvtx_chipID[daughter_id].push_back(chipId);
        intt_ladderZID[daughter_id].push_back(ladderZId);
        intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
        tpc_sectorID[daughter_id].push_back(sectorId);
        tpc_side[daughter_id].push_back(side);
      }
    }
  }

  if (tpcseed)
  {
    for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));

      unsigned int staveId;
      unsigned int chipId;
      unsigned int ladderZId;
      unsigned int ladderPhiId;
      unsigned int sectorId;
      unsigned int side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();

      if (Use["TPC"] && trackerID == TrkrDefs::tpcId)
      {
        sectorId = TpcDefs::getSectorId(cluster_key);
        side = TpcDefs::getSide(cluster_key);
        ++detector_nHits_TPC[daughter_id];
      }
      else if (Use["TPOT"] && trackerID == TrkrDefs::micromegasId)
      {
        ++detector_nHits_TPOT[daughter_id];
      }

      if (m_config.detailedTracking)
      {
        mvtx_staveID[daughter_id].push_back(staveId);
        mvtx_chipID[daughter_id].push_back(chipId);
        intt_ladderZID[daughter_id].push_back(ladderZId);
        intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
        tpc_sectorID[daughter_id].push_back(sectorId);
        tpc_side[daughter_id].push_back(side);
      }
    }
  }

  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter)
  {
    SvtxTrackState *tstate = state_iter->second;
    if (tstate->get_pathlength() != 0)  // The first track state is an extrapolation so has no cluster
    {
      auto stateckey = tstate->get_cluskey();
      TrkrCluster *cluster = dst_clustermap->findCluster(stateckey);
      if (!cluster)
      {
        // do not have associated cluster, could be track states projected to calo system
        continue;
      }
      auto global = geometry->getGlobalPosition(stateckey, cluster);

      if (m_config.detailedTracking)
      {
        residual_x[daughter_id].push_back(global.x() - tstate->get_x());
        residual_y[daughter_id].push_back(global.y() - tstate->get_y());
        residual_z[daughter_id].push_back(global.z() - tstate->get_z());
      }

      uint8_t id = TrkrDefs::getTrkrId(stateckey);

      switch (id)
      {
      case TrkrDefs::mvtxId:
        ++detector_nStates_MVTX[daughter_id];
        break;
      case TrkrDefs::inttId:
        ++detector_nStates_INTT[daughter_id];
        break;
      case TrkrDefs::tpcId:
        ++detector_nStates_TPC[daughter_id];
        break;
      case TrkrDefs::micromegasId:
        ++detector_nStates_TPOT[daughter_id];
        break;
      default:
        // std::cout << "Cluster key doesnt match a tracking system, could be related with projected track state to calorimeter system" << std::endl;
        break;
      }
    }
  }
}

void KFParticle_detectorInfo::clear(int nTracks)
{
  for (int i = 0; i < nTracks; ++i)
  {
    // Detector vectors
    residual_x[i].clear();
    residual_y[i].clear();
    residual_z[i].clear();
    detector_layer[i].clear();
    mvtx_staveID[i].clear();
    mvtx_chipID[i].clear();
    intt_ladderZID[i].clear();
    intt_ladderPhiID[i].clear();
    tpc_sectorID[i].clear();
    tpc_side[i].clear();

    detector_nStates_MVTX[i] = 0;
    detector_nStates_INTT[i] = 0;
    detector_nStates_TPC[i] = 0;
    detector_nStates_TPOT[i] = 0;
  }
}
