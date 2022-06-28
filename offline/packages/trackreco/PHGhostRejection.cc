#include "PHGhostRejection.h"

#include "PHGhostRejection.h"   

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase_historic/TrackSeed.h>     
#include <trackbase_historic/TrackSeedContainer.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair

//____________________________________________________________________________..
PHGhostRejection::PHGhostRejection(unsigned int verbosity)
  : m_verbosity(verbosity)
{
}

//____________________________________________________________________________..
PHGhostRejection::~PHGhostRejection()
{
}

//____________________________________________________________________________..
void PHGhostRejection::rejectGhostTracks(std::vector<float> &trackChi2)
{

  if(!m_trackMap || m_positions.size() == 0)
    {
      std::cout << "Missing containers, will not run TPC seed ghost rejection"
		<< std::endl;
      return;
    }

  if(m_verbosity > 0)
    std::cout << "PHGhostRejection beginning track map size " << m_trackMap->size() << std::endl;

  // Try to eliminate repeated tracks

  std::set<unsigned int> matches_set;
  std::multimap<unsigned int, unsigned int>  matches;
  
  for (unsigned int trid1 = 0;
       trid1 != m_trackMap->size(); 
       ++trid1)
    {
      TrackSeed* track1 = m_trackMap->get(trid1);
      if(!track1) 
	{ continue; }
      float track1phi = track1->get_phi(m_positions); 
      float track1x = track1->get_x();
      float track1y = track1->get_y();
      float track1z = track1->get_z();
      float track1eta = track1->get_eta();
      for (unsigned int trid2 = trid1;
	   trid2 != m_trackMap->size(); 
	   ++trid2)
	{
	  if(trid1  ==  trid2) continue;
	  
	  TrackSeed* track2 = m_trackMap->get(trid2);
	  if(!track2)
	    { continue; }
	  if(fabs( track1phi - track2->get_phi(m_positions)) < _phi_cut &&
	     fabs( track1eta - track2->get_eta() ) < _eta_cut &&
	     fabs( track1x - track2->get_x() ) < _x_cut &&
	     fabs( track1y - track2->get_y() ) < _y_cut &&
	     fabs( track1z - track2->get_z() ) < _z_cut
	      )
	    {
	      matches_set.insert(trid1);
	      matches.insert( std::pair( trid1, trid2) );

	      if(m_verbosity > 1)
		std::cout << "Found match for tracks " << trid1 << " and " << trid2 << std::endl;
	    }
	}
    }

  std::set<unsigned int> ghost_reject_list;

  for(auto set_it : matches_set)
    {
      if(ghost_reject_list.find(set_it) != ghost_reject_list.end()) continue;  // already rejected  

      auto match_list = matches.equal_range(set_it);

      auto tr1 = m_trackMap->get(set_it);
      double best_qual = trackChi2.at(set_it);
      unsigned int best_track = set_it;

      if(m_verbosity > 1)  
	std::cout << " ****** start checking track " << set_it << " with best quality " << best_qual << " best_track " << best_track << std::endl;

      for (auto it=match_list.first; it!=match_list.second; ++it)
	{
	  if(m_verbosity > 1)
	    std::cout << "    match of track " << it->first << " to track " << it->second << std::endl;
	  
	  auto tr2 = m_trackMap->get(it->second);	  

	  // Check that these two tracks actually share the same clusters, if not skip this pair

	  bool is_same_track = checkClusterSharing(tr1, set_it,
						   tr2, it->second);
	  if(!is_same_track) continue;

	  // which one has the best quality?
	  double tr2_qual = trackChi2.at(it->second);
	  if(m_verbosity > 1)
	    {
	      std::cout << "       Compare: best quality " << best_qual << " track 2 quality " << tr2_qual << std::endl;
	      std::cout << "       tr1: phi " << tr1->get_phi(m_positions) << " eta " << tr1->get_eta() 
			<<  " x " << tr1->get_x() << " y " << tr1->get_y() << " z " << tr1->get_z() << std::endl;
	      std::cout << "       tr2: phi " << tr2->get_phi(m_positions) << " eta " << tr2->get_eta() 
			<<  " x " << tr2->get_x() << " y " << tr2->get_y() << " z " << tr2->get_z() << std::endl;
	    }

	  if(tr2_qual < best_qual)
	    {
	      if(m_verbosity > 1)
		std::cout << "       --------- Track " << it->second << " has better quality, erase track " << best_track << std::endl;
	      ghost_reject_list.insert(best_track);
	      best_qual = tr2_qual;
	      best_track = it->second;
	    }
	  else
	    {
	      if(m_verbosity > 1)
		std::cout << "       --------- Track " << best_track << " has better quality, erase track " << it->second << std::endl;
	      ghost_reject_list.insert(it->second);
	    }

	}
      if(m_verbosity > 1)
	std::cout << " best track " << best_track << " best_qual " << best_qual << std::endl;      

    }

  // delete ghost tracks
  for(auto it : ghost_reject_list)
    {
      if(m_verbosity > 1)
	std::cout << " erasing track ID " << it << std::endl;

      m_trackMap->erase(it);
    }
  
  if(m_verbosity > 0)
    std::cout << "Track map size after deleting ghost tracks: " << m_trackMap->size() << std::endl;


  return;
}

bool PHGhostRejection::checkClusterSharing(TrackSeed* tr1, unsigned int trid1,
					   TrackSeed* tr2, unsigned int trid2)
{
  // Check that tr1 and tr2 share many clusters

  bool is_same_track = false;

  std::multimap<TrkrDefs::cluskey, unsigned int> cluskey_map;
  std::vector<TrkrDefs::cluskey> clusterkeys;

  for (TrackSeed::ConstClusterKeyIter key_iter = tr1->begin_cluster_keys();
       key_iter != tr1->end_cluster_keys();
       ++key_iter)
    {
	  TrkrDefs::cluskey cluster_key = *key_iter;
	  if(m_verbosity > 2) std::cout << " track id: " << trid1 <<  " adding clusterkey " << cluster_key << std::endl;
	  cluskey_map.insert( std::make_pair(cluster_key, trid1) );
	  clusterkeys.push_back(cluster_key);
    }

  for (TrackSeed::ConstClusterKeyIter key_iter = tr2->begin_cluster_keys();
       key_iter != tr2->end_cluster_keys();
       ++key_iter)
    {
	  TrkrDefs::cluskey cluster_key = *key_iter;  
	  if(m_verbosity > 2) std::cout << " track id: " << trid2 <<  " adding clusterkey " << cluster_key << std::endl;
	  cluskey_map.insert( std::make_pair(cluster_key, trid2) );
    }
  
  unsigned int nclus = clusterkeys.size();

  unsigned int nclus_used = 0;
  for (TrkrDefs::cluskey& cluskey : clusterkeys)
  {
    if(cluskey_map.count(cluskey)>0) nclus_used++;
  }

  if( (float) nclus_used / (float) nclus > 0.5)
    is_same_track = true;

  if(m_verbosity > 1)
    std::cout << " tr1 " << trid1 << " tr2 " << trid2 << " nclus_used " << nclus_used << " nclus " << nclus << std::endl;
  if(m_verbosity > 0)
    if(!is_same_track) std::cout << "   ***** not the same track! ********" << " tr1 " << trid1 << " tr2 " 
				 << trid2 << " nclus_used " << nclus_used << " nclus " << nclus << std::endl;

  return is_same_track;
}
