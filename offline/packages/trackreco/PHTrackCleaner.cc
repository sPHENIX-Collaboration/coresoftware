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

  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size() 
	      << " _seed_track_map size " << _seed_track_map->size() << std::endl;

  std::set<unsigned int> track_keep_list;
  std::set<unsigned int> track_delete_list;

  unsigned int good_track = 0;  // for diagnostic output only
  unsigned int ok_track = 0;   // tracks to keep

  // loop over the TPC seed - track map and make a set containing all TPC seed ID's
  std::set<unsigned int> seed_id_list;
  auto map_range =  _seed_track_map->getAll();
  for(auto it = map_range.first; it != map_range.second; ++it)
    {
      seed_id_list.insert( (*it).first );
    }

  if(Verbosity() > 0)
    std::cout << " seed_id_list size " << seed_id_list.size() << std::endl;

  // loop over the TPC seed ID's

  for(auto seed_iter = seed_id_list.begin(); seed_iter != seed_id_list.end(); ++seed_iter)
    {
      unsigned int tpc_id = *seed_iter;

      if(Verbosity() > 1)
	std::cout << " TPC ID " << tpc_id << std::endl;

      auto tpc_range =   _seed_track_map->getAssocTracks(tpc_id);

      unsigned int best_id = 99999;
      double min_chisq = 99999.0;
      unsigned int best_ndf = 1;
      for (auto it = tpc_range.first; it !=tpc_range.second; ++it)
	{
	  unsigned int track_id = it->second;
	  
	  // note that the track no longer exists if it failed in the Acts fitter
	  _track = _track_map->get(track_id);
	  if(_track)
	    {
	      if(Verbosity() > 1)	      
		std::cout << "        track ID " << track_id << " chisq " << _track->get_chisq() << " ndf " << _track->get_ndf() << " min_chisq " << min_chisq << std::endl;

	      // only accept tracks with nclus > min_clusters
	      if(_track->get_chisq() < min_chisq && _track->size_cluster_keys() > min_clusters)
		{
		  min_chisq = _track->get_chisq();
		  best_id = track_id;
		  best_ndf = _track->get_ndf();
		}
	    }
	}

      if(best_id != 99999)
	{
	  double qual = min_chisq / best_ndf;

	  if(Verbosity() > 1)
	    std::cout << "        best track for tpc_id " << tpc_id << " has track_id " << best_id << " chisq " << min_chisq << " chisq/ndf " << qual << std::endl;

	  if(qual < 30)
	    {
	      track_keep_list.insert(best_id);
	      ok_track++;
	      if(qual < 10.0)
		good_track++;
	    }
	}
      else
	{
	  if(Verbosity() > 1)
	    std::cout << "        no track exists  for tpc_id " << tpc_id << std::endl;
	}
    }

  if(Verbosity() > 0)
    std::cout << " Number of good tracks with qual < 10  is " << good_track << " OK tracks " << ok_track << std::endl; 

  // make a list of tracks that did not make the keep list
  for(auto track_it = _track_map->begin(); track_it != _track_map->end(); ++track_it)
    {
      auto id = track_it->first;

      auto set_it = track_keep_list.find(id);
      if(set_it == track_keep_list.end())
	{
	  if(Verbosity() > 1)
	    std::cout << "    add id " << id << " to track_delete_list " << std::endl;
	  track_delete_list.insert(id);
	}
    }

  if(Verbosity() > 0)
    std::cout << " track_delete_list size " << track_delete_list.size() << std::endl;

  // delete failed tracks
  for(auto it = track_delete_list.begin(); it != track_delete_list.end(); ++it)
    {
      if(Verbosity() > 1)
	std::cout << " erasing track ID " << *it << std::endl;
      _track_map->erase(*it);
    }

  if(Verbosity() > 0)
    std::cout << "Track map size after choosing best silicon match: " << _track_map->size() << std::endl;

  // now we have a single silicon match per TPC seed
  // Try to eliminate tracks based on repeated TPC seeds

  if(_reject_ghosts)
    {
      findGhostTracks();

      if(Verbosity() > 0)
	std::cout << "Track map size after deleting ghost tracks: " << _track_map->size() << std::endl;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackCleaner::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHTrackCleaner::GetNodes(PHCompositeNode* topNode)
{

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _seed_track_map = findNode::getClass<TpcSeedTrackMap>(topNode, "TpcSeedTrackMap");
  if (!_seed_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TpcSeedTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHTrackCleaner::findGhostTracks()
{
  std::set<unsigned int> matches_set;
  std::multimap<unsigned int, unsigned int>  matches;
  for (auto tr1_iter = _track_map->begin();
       tr1_iter != _track_map->end(); 
       ++tr1_iter)
    {
      auto track1 = (tr1_iter)->second;
      
      for (auto tr2_iter = tr1_iter;
	   tr2_iter != _track_map->end(); 
	   ++tr2_iter)
	{
	  if((tr2_iter)->first  ==  (tr1_iter)->first) continue;
	  
	  auto track2 = (tr2_iter)->second;
	  if( 
	     fabs( track1->get_phi() - track2->get_phi() ) < _phi_cut &&
	     fabs( track1->get_eta() - track2->get_eta() ) < _eta_cut &&
	     fabs( track1->get_x() - track2->get_x() ) < _x_cut &&
	     fabs( track1->get_y() - track2->get_y() ) < _y_cut &&
	     fabs( track1->get_z() - track2->get_z() ) < _z_cut
	      )
	    {
	      matches_set.insert(tr1_iter->first);
	      matches.insert( std::pair( (tr1_iter)->first, (tr2_iter)->first) );

	      if(Verbosity() > 0)
		std::cout << "Found match for tracks " << (tr1_iter)->first << " and " << (tr2_iter)->first << std::endl;
	    }
	}
    }

  std::set<unsigned int> ghost_reject_list;

  for(auto set_it : matches_set)
    {
      auto match_list = matches.equal_range(set_it);

      auto tr1 = _track_map->get(set_it);
      double best_qual = tr1->get_chisq() / tr1->get_ndf();
      unsigned int best_track = set_it;

      if(Verbosity() > 0)  
	std::cout << " start checking track " << set_it << " with best quality " << best_qual << " best_track " << best_track << std::endl;

      for (auto it=match_list.first; it!=match_list.second; ++it)
	{
	  if(Verbosity() > 0)
	    std::cout << "    match of track " << it->first << " to track " << it->second << std::endl;
	  
	  // which one has best quality?
	  auto tr2 = _track_map->get(it->second);	  
	  double tr2_qual = tr2->get_chisq() / tr2->get_ndf();
	  if(Verbosity() > 0)
	    std::cout << "       best quality " << best_qual << " track 2 quality " << tr2_qual << std::endl;

	  if(tr2_qual < best_qual)
	    {
	      ghost_reject_list.insert(best_track);
	      best_qual = tr2_qual;
	      best_track = it->second;
	    }
	  else
	    ghost_reject_list.insert(it->second);

	}
      if(Verbosity() > 0)
	std::cout << " best track " << best_track << " best_qual " << best_qual << std::endl;      

    }

  // delete ghost tracks
  //for(auto it = ghost_reject_list.begin(); it != ghost_reject_list.end(); ++it)
  for(auto it : ghost_reject_list)
    {
      if(Verbosity() > 1)
	std::cout << " erasing track ID " << it << std::endl;

      _track_map->erase(it);
    }
  
  return;
}

