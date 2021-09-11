#include "PHGhostRejection.h"

#include "PHGhostRejection.h"   

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
PHGhostRejection::PHGhostRejection(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
PHGhostRejection::~PHGhostRejection()
{

}

//____________________________________________________________________________..
int PHGhostRejection::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHGhostRejection::process_event(PHCompositeNode * /*topNode*/)
{

  if(Verbosity() > 0)
    std::cout << PHWHERE << " Beginning track map size " << _track_map->size() << std::endl;

  // Try to eliminate repeated tracks

  findGhostTracks();

  if(Verbosity() > 0)
    std::cout << "Track map size after deleting ghost tracks: " << _track_map->size() << std::endl;

  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHGhostRejection::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHGhostRejection::GetNodes(PHCompositeNode* topNode)
{

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHGhostRejection::findGhostTracks()
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

	      if(Verbosity() > 1)
		std::cout << "Found match for tracks " << (tr1_iter)->first << " and " << (tr2_iter)->first << std::endl;
	    }
	}
    }

  std::set<unsigned int> ghost_reject_list;

  for(auto set_it : matches_set)
    {
      if(ghost_reject_list.find(set_it) != ghost_reject_list.end()) continue;  // already rejected  

      auto match_list = matches.equal_range(set_it);

      auto tr1 = _track_map->get(set_it);
      double best_qual = tr1->get_chisq() / tr1->get_ndf();
      unsigned int best_track = set_it;

      if(Verbosity() > 1)  
	std::cout << " ****** start checking track " << set_it << " with best quality " << best_qual << " best_track " << best_track << std::endl;

      for (auto it=match_list.first; it!=match_list.second; ++it)
	{
	  if(Verbosity() > 1)
	    std::cout << "    match of track " << it->first << " to track " << it->second << std::endl;
	  
	  // which one has best quality?
	  auto tr2 = _track_map->get(it->second);	  
	  double tr2_qual = tr2->get_chisq() / tr2->get_ndf();
	  if(Verbosity() > 1)
	    {
	      std::cout << "       Compare: best quality " << best_qual << " track 2 quality " << tr2_qual << std::endl;
	      std::cout << "       tr1: phi " << tr1->get_phi() << " eta " << tr1->get_eta() 
			<<  " x " << tr1->get_x() << " y " << tr1->get_y() << " z " << tr1->get_z() << std::endl;
	      std::cout << "       tr2: phi " << tr2->get_phi() << " eta " << tr2->get_eta() 
			<<  " x " << tr2->get_x() << " y " << tr2->get_y() << " z " << tr2->get_z() << std::endl;
	    }

	  if(tr2_qual < best_qual)
	    {
	      if(Verbosity() > 1)
		std::cout << "       --------- Track " << it->second << " has better quality, erase track " << best_track << std::endl;
	      ghost_reject_list.insert(best_track);
	      best_qual = tr2_qual;
	      best_track = it->second;
	    }
	  else
	    {
	      if(Verbosity() > 1)
		std::cout << "       --------- Track " << best_track << " has better quality, erase track " << it->second << std::endl;
	      ghost_reject_list.insert(it->second);
	    }

	}
      if(Verbosity() > 1)
	std::cout << " best track " << best_track << " best_qual " << best_qual << std::endl;      

    }

  // delete ghost tracks
  for(auto it : ghost_reject_list)
    {
      if(Verbosity() > 1)
	std::cout << " erasing track ID " << it << std::endl;

      _track_map->erase(it);
    }
  
  return;
}

