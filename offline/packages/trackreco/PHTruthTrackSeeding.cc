#include "PHTruthTrackSeeding.h"

#include "AssocInfoContainer.h"

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTra...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, Svtx...
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl
#include <map>       // for multimap, map<>::c...
#include <memory>
#include <utility>  // for pair

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

class PHCompositeNode;

using namespace std;

PHTruthTrackSeeding::PHTruthTrackSeeding(const std::string& name)
  : PHTrackSeeding(name)
  , _g4truth_container(nullptr)
  , phg4hits_tpc(nullptr)
  , phg4hits_intt(nullptr)
  , phg4hits_mvtx(nullptr)
  , hittruthassoc(nullptr)
  , clusterhitassoc(nullptr)
  , _seeding_layers({7, 13, 19, 25, 31, 37, 40})
  , _min_clusters_per_track(0)
{
}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode)
{
  cout << "Enter PHTruthTrackSeeding:: Setup" << endl;

  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process()
{
  typedef std::map<int, std::set<TrkrCluster*> > TrkClustersMap;
  TrkClustersMap m_trackID_clusters;

  // loop over all clusters
  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();
  for (TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
  {
    TrkrCluster* cluster = clusiter->second;
    TrkrDefs::cluskey cluskey = clusiter->first;
    unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);

    // get the hits for this cluster
    TrkrClusterHitAssoc::ConstRange hitrange = clusterhitassoc->getHits(cluskey);  // returns range of pairs {cluster key, hit key} for this cluskey
    for (TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
    {
      TrkrDefs::hitkey hitkey = clushititer->second;
      // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);

      // get all of the g4hits for this hitkey
      std::multimap<TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;
      hittruthassoc->getG4Hits(hitsetkey, hitkey, temp_map);  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
      for (std::multimap<TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter = temp_map.begin(); htiter != temp_map.end(); ++htiter)
      {
        // extract the g4 hit key here and add the hits to the set
        PHG4HitDefs::keytype g4hitkey = htiter->second.second;
        PHG4Hit* phg4hit;
        if (trkrid == TrkrDefs::tpcId)
          phg4hit = phg4hits_tpc->findHit(g4hitkey);
        else if (trkrid == TrkrDefs::inttId)
          phg4hit = phg4hits_intt->findHit(g4hitkey);
        else
          phg4hit = phg4hits_mvtx->findHit(g4hitkey);

        int particle_id = phg4hit->get_trkid();

        TrkClustersMap::iterator it = m_trackID_clusters.find(particle_id);

        if (it != m_trackID_clusters.end())
        {
          it->second.insert(cluster);
        }
        else
        {
          std::set<TrkrCluster*> clusters;
          clusters.insert(cluster);
          m_trackID_clusters.insert(std::pair<int, std::set<TrkrCluster*> >(particle_id, clusters));
        }
      }  // loop over g4hits associated with hit
    }    // loop over hits associated with cluster
  }      // loop over clusters

  //==================================

  // Build track
  for (TrkClustersMap::const_iterator trk_clusters_itr = m_trackID_clusters.begin();
       trk_clusters_itr != m_trackID_clusters.end(); ++trk_clusters_itr)
  {
    if (trk_clusters_itr->second.size() > _min_clusters_per_track)
    {
      std::unique_ptr<SvtxTrack_FastSim> svtx_track(new SvtxTrack_FastSim());

      svtx_track->set_id(_track_map->size());
      svtx_track->set_truth_track_id(trk_clusters_itr->first);

      // dummy values, set px to make it through the minimum pT cut
      svtx_track->set_px(10.);
      svtx_track->set_py(0.);
      svtx_track->set_pz(0.);
      for (TrkrCluster* cluster : trk_clusters_itr->second)
      {
        svtx_track->insert_cluster_key(cluster->getClusKey());
        _assoc_container->SetClusterTrackAssoc(cluster->getClusKey(), svtx_track->get_id());
      }
      _track_map->insert(svtx_track.get());
    }
  }

  if (Verbosity() >= 2)
  {
    cout << "Loop over SvtxTrackMap entries " << endl;
    for (SvtxTrackMap::Iter iter = _track_map->begin();
         iter != _track_map->end(); ++iter)
    {
      SvtxTrack* svtx_track = iter->second;

      //svtx_track->identify();
      //continue;

      cout << "Track ID: " << svtx_track->get_id() << ", Dummy Track pT: "
           << svtx_track->get_pt() << ", Truth Track/Particle ID: "
           << svtx_track->get_truth_track_id() << endl;

      //Print associated clusters;
      for (SvtxTrack::ConstClusterKeyIter iter =
               svtx_track->begin_cluster_keys();
           iter != svtx_track->end_cluster_keys(); ++iter)
      {
        TrkrDefs::cluskey cluster_key = *iter;
        TrkrCluster* cluster = _cluster_map->findCluster(cluster_key);
        float radius = sqrt(
            cluster->getX() * cluster->getX() + cluster->getY() * cluster->getY());
        cout << "       cluster ID: "
             << cluster->getClusKey() << ", cluster radius: " << radius
             << endl;
      }
    }
  }
  //==================================

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::GetNodes(PHCompositeNode* topNode)
{
  _g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_g4truth_container)
  {
    cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!clusterhitassoc)
  {
    cout << PHWHERE << "Failed to find TRKR_CLUSTERHITASSOC node, quit!" << endl;
    exit(1);
  }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    cout << PHWHERE << "Failed to find TRKR_HITTRUTHASSOC node, quit!" << endl;
    exit(1);
  }

  phg4hits_tpc = findNode::getClass<PHG4HitContainer>(
      topNode, "G4HIT_TPC");

  phg4hits_intt = findNode::getClass<PHG4HitContainer>(
      topNode, "G4HIT_INTT");

  phg4hits_mvtx = findNode::getClass<PHG4HitContainer>(
      topNode, "G4HIT_MVTX");

  if (!phg4hits_tpc and phg4hits_intt and !phg4hits_mvtx)
  {
    if (Verbosity() >= 0)
    {
      cerr << PHWHERE << " ERROR: No PHG4HitContainer found!" << endl;
    }
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::End()
{
  return 0;
}
