/*!
 * \file DSTClusterPruning.cc
 * \author Alex Patton <aopatton@mit.edu>
 */

#include "DSTClusterPruning.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterv5.h>
// include new cluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>
#include <utility>  //for pair

#include <TFile.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TTree.h>
//_____________________________________________________________________

//_____________________________________________________________________
DSTClusterPruning::DSTClusterPruning(const std::string& name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________
int DSTClusterPruning::Init(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "Writer Init end" << std::endl;
  }
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "DSTClusterPruning::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if (!evalNode)
  {
    // create
    std::cout << "DSTClusterPruning::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode("EVAL");
    dstNode->addNode(evalNode);
  }

  auto svtxNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svtxNode)
  {
    // create
    std::cout << "DSTTrackArrayReader::Init - SVTX node missing - creating" << std::endl;
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  auto trkrNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "TRKR"));
  if (!trkrNode)
  {
    // create
    std::cout << "DSTTrackArrayReader::Init - TRKR node missing - creating" << std::endl;
    trkrNode = new PHCompositeNode("TRKR");
    dstNode->addNode(trkrNode);
  }

  // make new cluster container
  auto clsNode = findNode::getClass<TrkrClusterContainer>(trkrNode, "TRKR_CLUSTER_SEED");
  if (!clsNode)
  {
    auto newClusterNode = new PHIODataNode<PHObject>(new TrkrClusterContainerv4, "TRKR_CLUSTER_SEED", "PHObject");
    trkrNode->addNode(newClusterNode);
  }

  auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK)
  {
    return res;
  }
  if (Verbosity() > 0)
  {
    std::cout << "Writer Init end" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTClusterPruning::process_event(PHCompositeNode* /*topNode*/)
{
  // make topNode run in Init
  // Init(topNode);
  //  load nodes
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << "::" << __func__ << "::" << __LINE__ << std::endl;
    std::cout << "DSTClusterPruning::process_event" << std::endl;
  }

  prune_clusters();
  if (Verbosity() > 3)
  {
    print_clusters();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int DSTClusterPruning::load_nodes(PHCompositeNode* topNode)
{
  // get necessary nodes
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // look for reduced cluster
  m_reduced_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");

  // m_reduced_track_map = findNode::getClass<SvtxTrackMap>(topNode, "ReducedTrackContainer");
  m_track_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  m_tpc_track_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  m_silicon_track_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void DSTClusterPruning::prune_clusters()
{
  // use this to create object that looks through both tracks and clusters and saves into new object
  // make sure clusters exist
  // make sure tracks exist
  if (!(m_cluster_map && m_reduced_cluster_map && m_track_seed_container && m_silicon_track_seed_container && m_tpc_track_seed_container))
  {
    return;
  }
  for (const auto& trackseed : *m_track_seed_container)
  {
    if (!trackseed)
    {
      continue;
    }

    unsigned int tpcIndex = trackseed->get_tpc_seed_index();
    unsigned int siliconIndex = trackseed->get_silicon_seed_index();

    TrackSeed* TPCSeed = nullptr;
    if (tpcIndex <= m_tpc_track_seed_container->end() - m_tpc_track_seed_container->begin())
    {
      TPCSeed = *(m_tpc_track_seed_container->begin() + tpcIndex);
    }
    TrackSeed* SiliconSeed = nullptr;
    if (siliconIndex <= m_silicon_track_seed_container->end() - m_silicon_track_seed_container->begin())
    {
      SiliconSeed = *(m_silicon_track_seed_container->begin() + siliconIndex);
    }

    if (TPCSeed)
    {
      if (Verbosity() > 1)
      {
        std::cout << "We are about to loop over cluster keys in TPC Seed" << std::endl;
        TPCSeed->identify();
      }
      for (auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        auto cluster = m_cluster_map->findCluster(cluster_key);
        if (!cluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
          continue;
        }
        if (!m_reduced_cluster_map->findCluster(cluster_key))
        {
          m_cluster = new TrkrClusterv5();
          m_cluster->CopyFrom(cluster);
          m_reduced_cluster_map->addClusterSpecifyKey(cluster_key, m_cluster);
        }
      }
    }

    if (SiliconSeed)
    {
      if (Verbosity() > 1)
      {
        std::cout << "We are about to loop over cluster keys in Silicon Seed" << std::endl;
        SiliconSeed->identify();
      }
      for (auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        auto cluster = m_cluster_map->findCluster(cluster_key);
        if (!cluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
          continue;
        }
        if (!m_reduced_cluster_map->findCluster(cluster_key))
        {
          m_cluster = new TrkrClusterv5();
          m_cluster->CopyFrom(cluster);
          m_reduced_cluster_map->addClusterSpecifyKey(cluster_key, m_cluster);
        }
      }
    }
  }
}

// fill original clusters
//_____________________________________________________________________
void DSTClusterPruning::fill_clusters()
{
  // use this to create object that looks through both tracks and clusters and saves into new object
  // make sure clusters exist
  // if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;
  // make sure tracks exist
  if (!(m_cluster_map && m_reduced_cluster_map && m_track_seed_container && m_silicon_track_seed_container && m_tpc_track_seed_container))
  {
    return;
  }

  for (const auto& trackseed : *m_track_seed_container)
  {
    if (!trackseed)
    {
      continue;
    }

    unsigned int tpcIndex = trackseed->get_tpc_seed_index();
    unsigned int siliconIndex = trackseed->get_silicon_seed_index();

    TrackSeed* TPCSeed = nullptr;
    if (tpcIndex <= m_tpc_track_seed_container->end() - m_tpc_track_seed_container->begin())
    {
      TPCSeed = *(m_tpc_track_seed_container->begin() + tpcIndex);
    }
    TrackSeed* SiliconSeed = nullptr;
    if (siliconIndex <= m_silicon_track_seed_container->end() - m_silicon_track_seed_container->begin())
    {
      SiliconSeed = *(m_silicon_track_seed_container->begin() + siliconIndex);
    }

    if (TPCSeed)
    {
      for (auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        auto cluster = m_reduced_cluster_map->findCluster(cluster_key);
        if (!cluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
          continue;
        }
        if (!m_cluster_map->findCluster(cluster_key))
        {
          m_cluster = new TrkrClusterv5();
          m_cluster->CopyFrom(cluster);
          m_cluster_map->addClusterSpecifyKey(cluster_key, m_cluster);
        }
      }
    }

    if (SiliconSeed)
    {
      // std::cout << "We are about to loop over cluster keys in Silicon Seed" << std::endl;
      // SiliconSeed->identify();
      for (auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        auto cluster = m_reduced_cluster_map->findCluster(cluster_key);
        if (!cluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
          continue;
        }
        if (!m_cluster_map->findCluster(cluster_key))
        {
          m_cluster = new TrkrClusterv5();
          m_cluster->CopyFrom(cluster);
          m_cluster_map->addClusterSpecifyKey(cluster_key, m_cluster);

          // m_reduced_cluster_map->addClusterSpecifyKey(cluster_key, cluster);
        }
      }
    }
  }
}

// print clusters to test
void DSTClusterPruning::print_clusters()
{
  // make sure tracks exist
  if (!(m_track_map && m_cluster_map && m_reduced_cluster_map))
  {
    return;
  }

  for (const auto& trackpair : *m_track_map)
  {
    std::cout << "start of loop" << "\n";

    // unsigned int key = trackpair.first;
    const auto track = trackpair.second;

    TrackSeed* TPCSeed = track->get_tpc_seed();
    TrackSeed* SiliconSeed = track->get_silicon_seed();

    if (TPCSeed)
    {
      std::cout << "We are about to loop over cluster keys in TPC Seed" << std::endl;
      TPCSeed->identify();
      for (auto key_iter = TPCSeed->begin_cluster_keys(); key_iter != TPCSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        auto cluster = m_cluster_map->findCluster(cluster_key);
        auto reducedcluster = m_reduced_cluster_map->findCluster(cluster_key);
        if (!cluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
          continue;
        }
        if (!reducedcluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find reducedcluster for key " << cluster_key << std::endl;
          continue;
        }
        std::cout << "ClusterKey: " << cluster_key << std::endl;
        std::cout << "Cluster map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
        std::cout << "Reduced map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
      }
    }

    if (SiliconSeed)
    {
      std::cout << "We are about to loop over cluster keys in Silicon Seed" << std::endl;
      SiliconSeed->identify();
      for (auto key_iter = SiliconSeed->begin_cluster_keys(); key_iter != SiliconSeed->end_cluster_keys(); ++key_iter)
      {
        const auto& cluster_key = *key_iter;
        auto cluster = m_cluster_map->findCluster(cluster_key);

        auto reducedcluster = m_reduced_cluster_map->findCluster(cluster_key);
        if (!cluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find cluster for key " << cluster_key << std::endl;
          continue;
        }
        if (!reducedcluster)
        {
          std::cout << "DSTClusterPruning::evaluate_tracks - unable to find reducedcluster for key " << cluster_key << std::endl;
          continue;
        }
        std::cout << "ClusterKey: " << cluster_key << std::endl;
        std::cout << "Cluster map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
        std::cout << "Reduced map Cluster Local X: " << cluster->getLocalX() << " Local Y: " << cluster->getLocalY() << std::endl;
      }
    }

    std::cout << "end of loop" << "\n";
  }
}
