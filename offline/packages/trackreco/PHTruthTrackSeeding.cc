#include "PHTruthTrackSeeding.h"

#include "AssocInfoContainer.h"

//#include <trackbase_historic/SvtxClusterMap.h>
//#include <trackbase_historic/SvtxHitMap.h>
//#include <trackbase_historic/SvtxHit_v1.h>

#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex_v1.h>

#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <memory>

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << std::endl

using namespace std;

PHTruthTrackSeeding::PHTruthTrackSeeding(const std::string& name)
  : PHTrackSeeding(name)
  , _g4truth_container(nullptr)
  , phg4hits_tpc(nullptr)
  , phg4hits_intt(nullptr)
  , phg4hits_mvtx(nullptr)
    //  , hitsmap(nullptr)
    //, cells_svtx(nullptr)
    //, cells_intt(nullptr)
    //, cells_maps(nullptr)
  , _seeding_layers({7, 13, 19, 25, 31, 37, 40})
  , _min_clusters_per_track(0)
{
}

int PHTruthTrackSeeding::Setup(PHCompositeNode* topNode)
{
  int ret = PHTrackSeeding::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTruthTrackSeeding::Process()
{
  //==================================
  // The following has to be completely replaced
  //   need to get all clusters
  //   get all hits for each clusrer
  //   get all g4hits for each hit
  //   make a set named clusters to capture all unique cluster ID's
  //   redefine TrkClustersMap as  typedef std::map<int, std::set<TrkrCluster*> > TrkClustersMap  (maybe best to use clutsrkey?)
  //   make a map named m_trackid of the form:  std::pair<int, std::set<TrkrCluster*> >(g4hit particle_id, clusters))

  //----------------------------
  // new

  typedef std::map<int, std::set<TrkrCluster*> > TrkClustersMap;
  TrkClustersMap m_trackID_clusters;
  
  // loop over all clusters
  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();
  for(TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
    {
      TrkrCluster *cluster = clusiter->second;
      TrkrDefs::cluskey cluskey = clusiter->first;
      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);

      // get the hits for this cluster
     TrkrClusterHitAssoc::ConstRange hitrange = clusterhitassoc->getHits(cluskey);  // returns range of pairs {cluster key, hit key} for this cluskey
      for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
	{
	  TrkrDefs::hitkey hitkey = clushititer->second;
	  // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
	  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);	  
	  //if(layer < 7) cout << "  hitkey " << hitkey << endl;

	  // get all of the g4hits for this hitkey
	  std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
	  hittruthassoc->getG4Hits(hitsetkey, hitkey, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
	  for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); htiter != temp_map.end(); ++htiter) 
	    {
	      // extract the g4 hit key here and add the hits to the set
	      PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	      PHG4Hit * phg4hit;
	      if(trkrid == TrkrDefs::tpcId)
		phg4hit = phg4hits_tpc->findHit(g4hitkey);
	      else if(trkrid == TrkrDefs::inttId)
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
	}  // loop over hits associated with cluster
    }  // loop over clusters


  //-----------------------------------------------
  // old
  // Build TrackID -> Clusters map
  // typedef std::map<int, std::set<SvtxCluster*> > TrkClustersMap;
  // TrkClustersMap m_trackID_clusters;

  // for (SvtxClusterMap::ConstIter cluster_itr = _cluster_map->begin();
  //      cluster_itr != _cluster_map->end(); ++cluster_itr)
  // {
  //   SvtxCluster* cluster = cluster_itr->second;

  //   if (_seeding_layers.size() > 0 and
  //       (_seeding_layers.find(cluster->get_layer()) == _seeding_layers.end()))
  //     continue;

  //   SvtxHit* svtxhit = hitsmap->find(*cluster->begin_hits())->second;
  //   PHG4Cell* cell = nullptr;

  //   if (cells_svtx) cell = cells_svtx->findCell(svtxhit->get_cellid());
  //   if (!cell and cells_intt) cell = cells_intt->findCell(svtxhit->get_cellid());
  //   if (!cell and cells_maps) cell = cells_maps->findCell(svtxhit->get_cellid());

  //   if (!cell)
  //   {
  //     if (Verbosity() >= 1)
  //     {
  //       LogError("!cell");
  //     }
  //     continue;
  //   }

  //   //cell->identify();

  //   for (PHG4Cell::EdepConstIterator hits_it = cell->get_g4hits().first;
  //        hits_it != cell->get_g4hits().second; hits_it++)
  //   {
  //     PHG4Hit* phg4hit = nullptr;
  //     if (phg4hits_svtx) phg4hit = phg4hits_svtx->findHit(hits_it->first);
  //     if (!phg4hit and phg4hits_intt) phg4hit = phg4hits_intt->findHit(hits_it->first);
  //     if (!phg4hit and phg4hits_maps) phg4hit = phg4hits_maps->findHit(hits_it->first);

  //     if (!phg4hit)
  //     {
  //       if (Verbosity() >= 1)
  //       {
  //         LogError("!phg4hit");
  //       }
  //       continue;
  //     }

  //     //phg4hit->identify();

  //     int particle_id = phg4hit->get_trkid();

  //     TrkClustersMap::iterator it = m_trackID_clusters.find(particle_id);

  //     if (it != m_trackID_clusters.end())
  //     {
  //       it->second.insert(cluster);
  //     }
  //     else
  //     {
  //       std::set<SvtxCluster*> clusters;
  //       clusters.insert(cluster);
  //       m_trackID_clusters.insert(std::pair<int, std::set<SvtxCluster*> >(particle_id, clusters));
  //     }
  //   }
  // }

  // //--------------------------
  // // end old
  // //==================================
  // 

  //==================================
  // for this part, just replace Svtx types with Trkr types?

  // Build track
  for (TrkClustersMap::const_iterator trk_clusers_itr = m_trackID_clusters.begin();
       trk_clusers_itr != m_trackID_clusters.end(); ++trk_clusers_itr)
  {
    if (trk_clusers_itr->second.size() > _min_clusters_per_track)
    {
      std::unique_ptr<SvtxTrack_FastSim> svtx_track(new SvtxTrack_FastSim());

      //TODO implement the track ID
      svtx_track->set_id(_track_map->size());
      svtx_track->set_truth_track_id(trk_clusers_itr->first);
      //to make through minimum pT cut
      svtx_track->set_px(10.);
      svtx_track->set_py(0.);
      svtx_track->set_pz(0.);
      for (TrkrCluster* cluster : trk_clusers_itr->second)
      {
        svtx_track->insert_cluster(cluster->getClusKey());
        _assoc_container->SetClusterTrackAssoc(cluster->getClusKey(), svtx_track->get_id());
      }
      _track_map->insert(svtx_track.get());
    }
  }

  if (Verbosity() >= 2)
  {
    for (SvtxTrackMap::Iter iter = _track_map->begin();
         iter != _track_map->end(); ++iter)
    {
      SvtxTrack* svtx_track = iter->second;
      svtx_track->identify();
      continue;
      //Print associated clusters;
      for (SvtxTrack::ConstClusterIter iter =
               svtx_track->begin_clusters();
           iter != svtx_track->end_clusters(); ++iter)
      {
        unsigned int cluster_id = *iter;
	TrkrCluster* cluster = _cluster_map->findCluster(cluster_id);
        float radius = sqrt(
            cluster->getPosition(0) * cluster->getPosition(0) + cluster->getPosition(1) * cluster->getPosition(1));
        cout << "Track ID: " << svtx_track->get_id() << ", Track pT: "
             << svtx_track->get_pt() << ", Particle ID: "
             << svtx_track->get_truth_track_id() << ", cluster ID: "
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

  clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode,"TRKR_CLUSTERHITASSOC");
  if(!clusterhitassoc) 
    {
      cout << PHWHERE << "Failed to find TRKR_CLUSTERHITASSOC node, quit!" << endl;
      exit(1);
    }
  
  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  if(!hittruthassoc) 
    {
      cout << PHWHERE << "Failed to find TRKR_HITTRUTHASSOC node, quit!" << endl;
      exit(1);
    }

  phg4hits_tpc = findNode::getClass<PHG4HitContainer>(
      topNode, "G4HIT_SVTX");

  phg4hits_intt = findNode::getClass<PHG4HitContainer>(
      topNode, "G4HIT_SILICON_TRACKER");

  phg4hits_mvtx = findNode::getClass<PHG4HitContainer>(
      topNode, "G4HIT_MAPS");

  if (!phg4hits_tpc and phg4hits_intt and !phg4hits_mvtx)
  {
    if (Verbosity() >= 0)
    {
      cerr << PHWHERE << " ERROR: No PHG4HitContainer found!" << endl;
    }
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /*
  hitsmap = nullptr;
  // get node containing the digitized hits
  hitsmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!hitsmap)
  {
    cout << PHWHERE << "ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  */
  /*
  cells_svtx = findNode::getClass<PHG4CellContainer>(
      topNode, "G4CELL_SVTX");

  cells_intt = findNode::getClass<PHG4CellContainer>(
      topNode, "G4CELL_SILICON_TRACKER");

  cells_maps = findNode::getClass<PHG4CellContainer>(
      topNode, "G4CELL_MAPS");

  if (!cells_svtx and !cells_intt and !cells_maps)
  {
    if (Verbosity() >= 0)
    {
      cerr << PHWHERE << " ERROR: No PHG4CellContainer found!" << endl;
    }
    return Fun4AllReturnCodes::ABORTRUN;
  }
  */

  return Fun4AllReturnCodes::EVENT_OK;
}
