#include "PHTrackCleaner.h"

#include "PHTrackCleaner.h"   

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TpcSeedTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

//____________________________________________________________________________..
PHTrackCleaner::PHTrackCleaner(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
PHTrackCleaner::~PHTrackCleaner()
{

}

//____________________________________________________________________________..
int PHTrackCleaner::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHTrackCleaner::process_event(PHCompositeNode *topNode)
{

  //if(Verbosity() > 0)
  std::cout << PHWHERE << " track map size " << _track_map->size() << " _seed_track_map size " << _seed_track_map_class->SeedTrackMap.size() << std::endl;

  // loop over the TPC seed - track map and make a set containing all TPC seed ID's
  std::set<unsigned int> seed_id_list;
  std::multimap<unsigned int, unsigned int>::iterator it;
  for(it = _seed_track_map_class->SeedTrackMap.begin(); it != _seed_track_map_class->SeedTrackMap.end(); ++it)
    {
      std::cout << " TPC seed ID " << (*it).first << " track ID " << (*it).second << std::endl;
      seed_id_list.insert( (*it).first );
    }

  std::cout << " seed_id_list size " << seed_id_list.size() << std::endl;


  /*
  // loop over the tracks
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      _track = phtrk_iter->second;
      
      if (Verbosity() >= 1)
	{
	  std::cout << std::endl
	    << __LINE__
	    << ": Processing track itrack: " << phtrk_iter->first
	    << ": nhits: " << _track-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _track->get_phi()
		    << std::endl;
	}

      // Get the TPC clusters for this track
      std::map<unsigned int, TrkrCluster*> tpc_clusters;
      std::vector<TrkrCluster*> clusters;

      for (SvtxTrack::ConstClusterKeyIter key_iter = _track->begin_cluster_keys();
	   key_iter != _track->end_cluster_keys();
	   ++key_iter)
	{
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  unsigned int trkrId = TrkrDefs::getTrkrId(cluster_key);
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  if(trkrId != TrkrDefs::tpcId) continue;  // we want only TPC clusters

	  // get the cluster
	  TrkrCluster *tpc_clus =  _cluster_map->findCluster(cluster_key);

	  tpc_clusters.insert(std::make_pair(layer, tpc_clus));
	  clusters.push_back(tpc_clus);

	  if(Verbosity() > 10) 
	    std::cout << "  TPC cluster in layer " << layer << " with position " << tpc_clus->getX() 
		      << "  " << tpc_clus->getY() << "  " << tpc_clus->getZ() << " outer_clusters.size() " << tpc_clusters.size() << std::endl;
	}

      // need at least 3 clusters to fit a circle
      if(tpc_clusters.size() < 3)
	{
	  if(Verbosity() > 3) std::cout << PHWHERE << "  -- skip this tpc track, not enough clusters " << std::endl; 
	  continue;  // skip to the next TPC track
	}

      // fit a circle to the clusters
      double R = 0;
      double X0 = 0;
      double Y0 = 0;
      CircleFitByTaubin(clusters, R, X0, Y0);
      if(Verbosity() > 10) 
	std::cout << " Fitted circle has R " << R << " X0 " << X0 << " Y0 " << Y0 << std::endl;

      // toss tracks for which the fitted circle could not have come from the vertex
      if(R < 30.0) continue;

      // get the straight line representing the z trajectory in the form of z vs radius
      double A = 0; double B = 0;
      line_fit(clusters, A, B);
      if(Verbosity() > 10) 
	std::cout << " Fitted line has A " << A << " B " << B << std::endl;

      // Now we need to move each cluster associated with this track to the readout layer radius
      for (auto clus_iter = tpc_clusters.begin();
	   clus_iter != tpc_clusters.end(); 
	   ++clus_iter)
	{
	  unsigned int layer = clus_iter->first;
	  TrkrCluster *cluster = clus_iter->second;
	 
	  // get circle position at target surface radius 
	  double target_radius = layer_radius[layer-7];
	  int ret = get_circle_circle_intersection(target_radius, R, X0, Y0, cluster->getX(), cluster->getY(), _x_proj, _y_proj);
	  if(ret == Fun4AllReturnCodes::ABORTEVENT) continue;  // skip to next cluster
	  // z projection is unique
	  _z_proj = B + A * target_radius;
	  
	  // get circle position at cluster radius	  
	  double cluster_radius = sqrt(cluster->getX() * cluster->getX() + cluster->getY() * cluster->getY());
	  ret = get_circle_circle_intersection(cluster_radius, R, X0, Y0, cluster->getX(), cluster->getY(), _x_start, _y_start);
	  if(ret == Fun4AllReturnCodes::ABORTEVENT) continue;  // skip to next cluster
	  // z projection is unique
	  _z_start = B + A * cluster_radius;
	  
	  // calculate dx, dy, dz along circle trajectory from cluster radius to surface radius
	  double xnew = cluster->getX() - (_x_start - _x_proj);
	  double ynew = cluster->getY() - (_y_start - _y_proj);
	  double znew = cluster->getZ() - (_z_start - _z_proj);
	  
	  // now move the cluster to the surface radius
	  cluster->setX(xnew);
	  cluster->setY(ynew);
	  cluster->setZ(znew);

	  if(Verbosity() > 0)
	    {
	      std::cout << "*** cluster_radius " << cluster_radius << " cluster x,y,z: " << cluster->getX() << "  " << cluster->getY() << "  " << cluster->getZ() << std::endl;
	      std::cout << "    projection_radius " << target_radius << " proj x,y,z: " << _x_proj << "  " << _y_proj << "  " << _z_proj << std::endl; 
	      std::cout << "    traj_start_radius " << cluster_radius << " start x,y,z: "<< _x_start << "  " << _y_start << "  " << _z_start << std::endl; 
	      std::cout << "    moved_clus_radius " << target_radius << " final x,y,z: "<< xnew << "  " << ynew << "  " << znew << std::endl; 
	    }
	  
	}
    }
  */
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackCleaner::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHTrackCleaner::GetNodes(PHCompositeNode* topNode)
{
  /*
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _seed_track_map_class = findNode::getClass<TpcSeedTrackMap>(topNode, "TpcSeedTrackMap");
  if (!_seed_track_map_class)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TpcSeedTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

