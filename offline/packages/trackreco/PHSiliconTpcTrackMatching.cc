#include "PHSiliconTpcTrackMatching.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TpcDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrClusterv3.h>   
#include <trackbase/TrkrClusterContainer.h>   
#include <trackbase/TrkrClusterCrossingAssoc.h>   

#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TF1.h>

#include <climits>                            // for UINT_MAX
#include <iostream>                            // for operator<<, basic_ostream
#include <cmath>                              // for fabs, sqrt
#include <set>                                 // for _Rb_tree_const_iterator
#include <utility>                             // for pair
#include <memory>

using namespace std;

//____________________________________________________________________________..
PHSiliconTpcTrackMatching::PHSiliconTpcTrackMatching(const std::string &name):
  SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

//____________________________________________________________________________..
PHSiliconTpcTrackMatching::~PHSiliconTpcTrackMatching()
{

}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::InitRun(PHCompositeNode *topNode)
{
  UpdateParametersWithMacro();

  // put these in the output file
  cout << PHWHERE << " Search windows: phi " << _phi_search_win << " eta " 
       << _eta_search_win << " _pp_mode " << _pp_mode << " _use_intt_time " << _use_intt_time << endl;


   int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//_____________________________________________________________________
void PHSiliconTpcTrackMatching::SetDefaultParameters()
{
  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html

  return;
}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::process_event(PHCompositeNode*)
{
  // _track_map contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // _svtx_seed_map contains the combined silicon and tpc track seeds

  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size() << " Silicon track map size " << _track_map_silicon->size() << endl;

  if(_track_map->size() == 0)
    return Fun4AllReturnCodes::EVENT_OK;
 
  // loop over the silicon seeds and add the crossing to them
  for (unsigned int trackid = 0; trackid != _track_map_silicon->size(); ++trackid) 
    {
      _tracklet_si = _track_map_silicon->get(trackid);	  
      if(!_tracklet_si) { continue; }

      short int crossing= getCrossingIntt(_tracklet_si);
      _tracklet_si->set_crossing(crossing);

      if(Verbosity() > 8)
	std::cout << " silicon stub: " << trackid << " eta " << _tracklet_si->get_eta()  << " pt " << _tracklet_si->get_pt()  << " si z " << _tracklet_si->get_z() << " crossing " << crossing << std::endl; 

      if(Verbosity() > 1) cout << " Si track " << trackid << " crossing " << crossing << endl;
    }  
  
  // Find all matches of tpc and si tracklets in eta and phi, x and y
  //     If _pp_mode is not set, a match in z is also required - gives same behavior as old code
  std::multimap<unsigned int, unsigned int> tpc_matches;
  std::set<unsigned int> tpc_matched_set;
  std::set<unsigned int> tpc_unmatched_set;
  findEtaPhiMatches(tpc_matched_set, tpc_unmatched_set, tpc_matches);

  // Check that the crossing number is consistent with the tracklet z mismatch, discard the match otherwise
  // Enabling this required a change to truth seeding, so that it sets the TPC seed z0 to the line fit value, not the truth
  checkCrossingMatches(tpc_matches);
  
  // We have a complete list of all eta/phi matched tracks in the map "tpc_matches"
  // make the combined track seeds from tpc_matches
  for(auto [tpcid, si_id] : tpc_matches)
    {
      auto svtxseed = std::make_unique<SvtxTrackSeed_v1>();
      svtxseed->set_silicon_seed_index(si_id);
      svtxseed->set_tpc_seed_index(tpcid);
      _svtx_seed_map->insert(svtxseed.get());
      
      if(Verbosity() > 1) std::cout << "  combined seed id " << _svtx_seed_map->size()-1 << " si id " << si_id << " tpc id " << tpcid  << std::endl;
    }

  // Also make the unmatched TPC seeds into SvtxTrackSeeds
  for(auto tpcid : tpc_unmatched_set)
    {
      auto svtxseed = std::make_unique<SvtxTrackSeed_v1>();
      svtxseed->set_tpc_seed_index(tpcid);
      _svtx_seed_map->insert(svtxseed.get());

      if(Verbosity() > 1) std::cout << "  converted unmatched TPC seed id " << _svtx_seed_map->size()-1 << " tpc id " << tpcid << std::endl;
    }
  /*
  // Future development: use the z-mismatch between the silicon and TPC tracklets to assign the crossing in case INTT clusters are missing
  // this will reuse some of the commented out methods at the end of this file
  tagMatchCrossing(tpc_matches, crossing_matches, tpc_crossing_map);
  */

  if(Verbosity() > 0)  
    {
      std::cout << "final svtx seed map size " << _svtx_seed_map->size() << std::endl;
    }

  if(Verbosity() > 1)
    { 
      for(const auto& seed : *_svtx_seed_map)
	seed->identify();

      cout << "PHSiliconTpcTrackMatching::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
    }

  return Fun4AllReturnCodes::EVENT_OK;

 }
  
int PHSiliconTpcTrackMatching::End(PHCompositeNode* )
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHSiliconTpcTrackMatching::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

   _cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!_cluster_crossing_map)
  {
    cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << endl;
    //return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_silicon = findNode::getClass<TrackSeedContainer>(topNode, _silicon_track_map_name);
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<TrackSeedContainer>(topNode, _track_map_name);
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name.c_str() << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _svtx_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if(!_svtx_seed_map)
    {
      std::cout << "Creating node SvtxTrackSeedContainer" << std::endl;
        /// Get the DST Node
      PHNodeIterator iter(topNode);
      PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  
      /// Check that it is there
      if (!dstNode)
	{
	  std::cerr << "DST Node missing, quitting" << std::endl;
	  throw std::runtime_error("failed to find DST node in PHActsSourceLinks::createNodes");
	}

      /// Get the tracking subnode
      PHNodeIterator dstIter(dstNode);
      PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));
      
      /// Check that it is there
      if (!svtxNode)
	{
	  svtxNode = new PHCompositeNode("SVTX");
	  dstNode->addNode(svtxNode);
	}
      
      _svtx_seed_map = new TrackSeedContainer_v1();
      PHIODataNode<PHObject> *node
	= new PHIODataNode<PHObject>(_svtx_seed_map, "SvtxTrackSeedContainer","PHObject");
      svtxNode->addNode(node);

    }

    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    if (!_cluster_map)
      {
	std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
	return Fun4AllReturnCodes::ABORTEVENT;
      }

   _tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
} 

void PHSiliconTpcTrackMatching::findEtaPhiMatches(  
		        std::set<unsigned int> &tpc_matched_set,
                        std::set<unsigned int> &tpc_unmatched_set,
			std::multimap<unsigned int, unsigned int> &tpc_matches )
{
  // loop over the TPC track seeds
  for (unsigned int phtrk_iter = 0;
       phtrk_iter < _track_map->size();
       ++phtrk_iter)
    {
      _tracklet_tpc = _track_map->get(phtrk_iter);
      if(!_tracklet_tpc) 
	{ continue; }
      
      unsigned int tpcid = phtrk_iter;
      if (Verbosity() > 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << tpcid
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet_tpc->get_phi(_cluster_map,_tGeometry)
	    << endl;
	}

      double tpc_phi = _tracklet_tpc->get_phi(_cluster_map,_tGeometry);
      double tpc_eta = _tracklet_tpc->get_eta();
      double tpc_pt = _tracklet_tpc->get_pt();
      if(Verbosity() > 8)
	std::cout << " tpc stub: " << tpcid << " eta " << tpc_eta << " pt " << tpc_pt << " tpc z " << _tracklet_tpc->get_z() << std::endl; 

      // this factor will increase the window size at low pT
      // otherwise the matching efficiency drops off at low pT
      // it would be better if this was a smooth function
      double mag = 1.0;
      if(tpc_pt < 6.0) mag = 2;
      if(tpc_pt < 3.0)  mag = 4.0;
      if(tpc_pt < 1.5)  mag = 6.0;

      if(Verbosity() > 3)
	{
	  cout << "TPC tracklet:" << endl;
	  _tracklet_tpc->identify();
	}

      double tpc_x = _tracklet_tpc->get_x();
      double tpc_y = _tracklet_tpc->get_y();
      double tpc_z = _tracklet_tpc->get_z();

      bool matched = false;

      // Now search the silicon track list for a match in eta and phi
      for (unsigned int phtrk_iter_si = 0;
	   phtrk_iter_si < _track_map_silicon->size(); 
	   ++phtrk_iter_si)
	{
	  _tracklet_si = _track_map_silicon->get(phtrk_iter_si);	  
	  if(!_tracklet_si)
	    { continue; }

	  bool eta_match = false;
	  double si_eta = _tracklet_si->get_eta();
	  if(  fabs(tpc_eta - si_eta) < _eta_search_win * mag)  eta_match = true;
//	  if(!eta_match) continue;
	  unsigned int siid = phtrk_iter_si;
	  double si_x = _tracklet_si->get_x();
	  double si_y = _tracklet_si->get_y();
	  double si_z = _tracklet_si->get_z();
	  bool position_match = false;
	  if(_pp_mode)
	    {
	      if(
		 fabs(tpc_x - si_x) < _x_search_win * mag
		 && fabs(tpc_y - si_y) < _y_search_win * mag 
		 )
		position_match = true;
	    }
	  else
	    {
	      if(
		 fabs(tpc_x - si_x) < _x_search_win * mag
		 && fabs(tpc_y - si_y) < _y_search_win * mag 
		 && fabs(tpc_z - si_z) < _z_search_win * mag 
		 )
		position_match = true;	    
	    }
	  
//	  if(!position_match)
//	    { continue; }

	  bool phi_match = false;
	  double si_phi = _tracklet_si->get_phi(_cluster_map,_tGeometry);
	  if(  fabs(tpc_phi - si_phi)  < _phi_search_win * mag) phi_match = true;
	  if(  fabs( fabs(tpc_phi - si_phi)  - 2.0 * M_PI)  < _phi_search_win * mag ) phi_match = true;
//	  if(!phi_match) continue;
	  if(Verbosity() > 3)
	    {
	      cout << " testing for a match for TPC track " << tpcid << " with pT " << _tracklet_tpc->get_pt() 
		   << " and eta " << _tracklet_tpc->get_eta() << " with Si track " << siid << " with crossing " << _tracklet_si->get_crossing() << endl;	  
	      cout << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " <<   tpc_phi-si_phi << " phi search " << _phi_search_win*mag  << " tpc_eta " << tpc_eta 
		   << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << " eta search " << _eta_search_win*mag << endl;
	      std::cout << "      tpc x " << tpc_x << " si x " << si_x << " tpc y " << tpc_y << " si y " << si_y << " tpc_z " << tpc_z  << " si z " << si_z << std::endl;
	      std::cout << "      x search " << _x_search_win*mag << " y search " << _y_search_win*mag << " z search " << _z_search_win*mag  << std::endl;
	    }

	  if(eta_match && phi_match && position_match)
	    {
	      // got a match, add to the list
	      // These stubs are matched in eta, phi, x and y already
	      matched = true;
	      tpc_matches.insert(std::make_pair(tpcid, siid));
	      tpc_matched_set.insert(tpcid);

	      if(Verbosity() > 1)  
		{
		  cout << " found a match for TPC track " << tpcid << " with Si track " << siid << endl;
		  cout << "          tpc_phi " << tpc_phi << " si_phi " <<  si_phi << " phi_match " << phi_match 
		       << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " eta_match " << eta_match << endl;
		  std::cout << "      tpc x " << tpc_x << " si x " << si_x << " tpc y " << tpc_y << " si y " << si_y << " tpc_z " << tpc_z  << " si z " << si_z << std::endl;
		}

	      // temporary!
	      if(_test_windows)
		cout << " Try_silicon:  pt " << tpc_pt << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi-si_phi  
		     << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << " tpc_x " << tpc_x << " tpc_y " << tpc_y << " tpc_z " << tpc_z 
		     << " dx " << tpc_x - si_x << " dy " << tpc_y - si_y << " dz " << tpc_z - si_z  
		     << endl;
	    }
	}
      // if no match found, keep tpc seed for fitting
      if(!matched)
        {
          if(Verbosity() > 1) cout << "inserted unmatched tpc seed " << tpcid << endl;
          tpc_unmatched_set.insert(tpcid);
        }
    }
  
  return;
}

short int PHSiliconTpcTrackMatching::getCrossingIntt(TrackSeed *si_track)
{
  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  
  std::vector<short int> intt_crossings = getInttCrossings(si_track);      
  
  bool keep_it = true;
  short int crossing_keep = 0;
  if(intt_crossings.size() == 0) 
    {
      keep_it = false ;
    }
  else
    {
      crossing_keep = intt_crossings[0];
      for(unsigned int ic=1; ic<intt_crossings.size(); ++ic)
	{	  
	  if(intt_crossings[ic] != crossing_keep)
	    {
	      if(Verbosity() > 1)
		{
		  std::cout << " Warning: INTT crossings not all the same " << " crossing_keep " << crossing_keep << " new crossing " << intt_crossings[ic] << " keep the first one in the list" << std::endl;  
		}
	    }
	}
    }
  
  if(keep_it)
    {            
      return crossing_keep;
    }

  return SHRT_MAX;
}

std::vector<short int> PHSiliconTpcTrackMatching::getInttCrossings(TrackSeed *si_track)
{
  std::vector<short int> intt_crossings;

  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  // loop over associated clusters to get keys for silicon cluster
  for (TrackSeed::ConstClusterKeyIter iter = si_track->begin_cluster_keys();
       iter != si_track->end_cluster_keys();
       ++iter)
    {
      
      TrkrDefs::cluskey cluster_key = *iter;
      const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);

      if(Verbosity() > 1) 
	{
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);

	  if(trkrid == TrkrDefs::mvtxId)
	    {
	      TrkrCluster *cluster =  _cluster_map->findCluster(cluster_key);	
	      if( !cluster ) continue;	  

	      Acts::Vector3 global  = _tGeometry->getGlobalPosition(cluster_key, cluster);

	      std::cout << "Checking  si Track " << _track_map_silicon->find(si_track) << " cluster " << cluster_key 
			<< " in layer " << layer  << " position " << global(0) << "  " << global(1) << "  " << global(2) 
			<< " eta " << si_track->get_eta() << std::endl;
	    }
	  else
	    std::cout << "Checking  si Track " << _track_map_silicon->find(si_track) << " cluster " << cluster_key  		
		      << " in layer " << layer  << " with eta " << si_track->get_eta() << std::endl;
	}      

      if(trkrid == TrkrDefs::inttId)
	{
	  TrkrCluster *cluster =  _cluster_map->findCluster(cluster_key);	
	  if( !cluster ) continue;	  
	  
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);
	  
	  // get the bunch crossings for all hits in this cluster
	  auto crossings = _cluster_crossing_map->getCrossings(cluster_key);
	  for(auto iter = crossings.first; iter != crossings.second; ++iter)
	    {
	      if(Verbosity() > 1) 
		std::cout << "                si Track " << _track_map_silicon->find(si_track) << " cluster " << iter->first << " layer " << layer << " crossing " << iter->second  << std::endl;
	      intt_crossings.push_back(iter->second);
	    }
	}
    }
  
  return intt_crossings;
}

void PHSiliconTpcTrackMatching::checkCrossingMatches( std::multimap<unsigned int, unsigned int> &tpc_matches )
{
  // if the  crossing was assigned correctly, the (crossing corrected) track position should satisfy the Z matching cut
  // this is a rough check that this is the case

  float vdrift = _tGeometry->get_drift_velocity();

  std::multimap<unsigned int, unsigned int> bad_map;

  for(auto [tpcid, si_id] : tpc_matches)
    {
      TrackSeed *tpc_track = _track_map->get(tpcid);
      TrackSeed *si_track = _track_map_silicon->get(si_id);
      short int crossing = si_track->get_crossing();

      if(crossing == SHRT_MAX) 
	{
	  if(Verbosity() > 2) {
	    std::cout << " drop si_track " << si_id << " with eta " << si_track->get_eta() << " and z " << si_track->get_z() << " because crossing is undefined " << std::endl; 
	  }
	  continue;
	}

      float z_si = si_track->get_z();
      float z_tpc = tpc_track->get_z();
      float z_mismatch = z_tpc-z_si;

      float mag_crossing_z_mismatch = fabs(crossing) * crossing_period * vdrift;

      // We do not know the sign  of the z mismatch for a given crossing unless we know the drift direction in the TPC, use magnitude
      // could instead look up any TPC cluster key in the track to get side
      // z-mismatch can occasionally be up to 2 crossings due to TPC extrapolation precision
      if( fabs( fabs(z_mismatch) - mag_crossing_z_mismatch ) < 3.0)
	{ 
	  if(Verbosity() > 1)	  
	    std::cout << "  Success:  crossing " << crossing << " tpcid " << tpcid << " si id " << si_id 
		      << " tpc z " << z_tpc << " si z " << z_si << " z_mismatch " << z_mismatch 
		      << " mag_crossing_z_mismatch " << mag_crossing_z_mismatch << " drift velocity " << vdrift << std::endl;
	}
      else
	{
	  if(Verbosity() > 1)
	    std::cout << "  FAILURE:  crossing " << crossing << " tpcid " << tpcid << " si id " << si_id 
		      << " tpc z " << z_tpc << " si z " << z_si << " z_mismatch " << z_mismatch 
		      << " mag_crossing_z_mismatch " << mag_crossing_z_mismatch << std::endl;

	  bad_map.insert(std::make_pair(tpcid, si_id));
	}
    }

  // remove bad entries from tpc_matches
  for(auto [tpcid, si_id] : bad_map)
    {
      // Have to iterate over tpc_matches and examine each pair to find the one matching bad_map
      // this logic works because we call the equal range on vertex_map for every id_pair
      // so we only delete one entry per equal range call
      auto ret = tpc_matches.equal_range(tpcid);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  if(it->first == tpcid && it->second == si_id)
	    {
	      if(Verbosity() > 1) 
		std::cout << "                        erasing tpc_matches entry for tpcid " << tpcid << " si_id " << si_id << std::endl;
	      tpc_matches.erase(it);
	      break;  // the iterator is no longer valid
	    }
	}
    }	  

  return;
}	  


/*
	{
	  // This section is to correct the TPC z positions of tracks for all bunch crossings without relying on INTT time  
	  //==================================================================================
	  
	  // All tracks treated as if we do not know the bunch crossing
	  // The crossing estimate here is crude, for now
	  tagMatchCrossing(tpc_matches, crossing_set, crossing_matches, tpc_crossing_map);
	  
	  // Sort candidates by the silicon tracklet Z position by putting them in si_sorted map
	  //      -- captures crossing, tpc_id and si_id
	  std::multimap<double, std::pair<unsigned int, unsigned int>> si_sorted_map;
	  for( auto ncross : crossing_set)
	    {
	      if(Verbosity() > 1) std::cout << " ncross = " << ncross << std::endl;
	      
	      auto ret  = crossing_matches. equal_range(ncross);            
	      for(auto it = ret.first; it != ret.second; ++it)
		{
		 if(Verbosity() > 1) 
		    std::cout << "  crossing " << it->first << " tpc_id " << it->second.first << " si_id " << it->second.second 
			      << " pT " << _track_map->get(it->second.first)->get_pt() 
			      << " eta " << _track_map->get(it->second.first)->get_eta() 
			      << " tpc_z " << _track_map->get(it->second.first)->get_z() 
			      << " si_z " << _track_map_silicon->get(it->second.second)->get_z() 
			      << std::endl;
		  
		  double z_si = _track_map_silicon->get(it->second.second)->get_z();
		  si_sorted_map.insert(std::make_pair(z_si, std::make_pair(it->second.first, it->second.second)));
		}	
	    }


	  // make a list of silicon vertices, and a multimap with all associated tpc/si pairs 
	  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>> vertex_map;
	  std::vector<double> vertex_list;
	  getSiVertexList(si_sorted_map, vertex_list, vertex_map);
	  
	  // find the crossing number for every track from the mismatch with the associated silicon vertex z position
	  std::map<unsigned int, short int> vertex_crossings_map;
	  getCrossingNumber(vertex_list, vertex_map, vertex_crossings_map);
	  
	  // Remove tracks from vertex_map where the vertex crossing is badly inconsistent with the initial crossing estimate
	  cleanVertexMap( vertex_crossings_map, vertex_map, tpc_crossing_map );
	  
	  // add silicon clusters to the surviving tracks
	  addSiliconClusters(vertex_map);

	  // add the crossing to the combined track
	  addTrackBunchCrossing(vertex_crossings_map, vertex_map);	  
	  
	  //===================================================
	}
*/

/*
void PHSiliconTpcTrackMatching::addTrackBunchCrossing(
						   std::map<unsigned int, short int> &vertex_crossings_map,
						   std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map)
{

 for(auto [ivert, crossing] : vertex_crossings_map)
    {
      if(Verbosity() > 1) 
	std::cout << "      ivert " << ivert << " crossing " << crossing << std::endl;		  

      // loop over all tracks associated with this vertex
      // remember that any tracks not associated with one of the vertices are unusable
      auto  ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  unsigned int tpcid = it->second.first;
	  SvtxTrack *tpc_track = _track_map->get(tpcid);
	  tpc_track->set_crossing(crossing);

	  if(Verbosity() > 1) std::cout << PHWHERE << "  Add bunch crossing to track " << tpcid << " crossing " << crossing << std::endl;		  
	}
    }
 return;
}
*/

/*
void PHSiliconTpcTrackMatching::getSiVertexList(  
						std::multimap<double, std::pair<unsigned int, unsigned int>> &si_sorted_map,
						  std::vector<double> &vertex_list,
						std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map)
{
  if(si_sorted_map.size() == 0) return;

  // process si_sorted_map to get the vertex list, and the average silicon z value for each vertex
  double zkeep = 999;
  double zavge = 0.0;
  double zwt = 0.0;
  unsigned int nvert = 0;
  for(auto it = si_sorted_map.begin(); it != si_sorted_map.end(); ++it)
    {
      double z_si = it->first;
      std::pair<unsigned int, unsigned int> id_pair = it->second;
      if(Verbosity() > 1) std::cout << " z_si " << z_si << " tpc_id " << it->second.first << " si_id " << it->second.second << std::endl; 
     
      if( (zkeep == 999) || (fabs(z_si - zkeep) < _si_vertex_dzmax) )
	{
	  vertex_map.insert(std::make_pair(nvert, id_pair));
	  zavge+= z_si;
	  zwt++;
	  zkeep = z_si;
	}
      else
	{
	  zavge /= zwt;
	  vertex_list.push_back(zavge);

	  // start a new vertex
	  nvert++;
	  zavge = z_si;
	  zwt = 1.0;
	  zkeep = z_si;
	  vertex_map.insert(std::make_pair(nvert, id_pair));	    
	}
    }
  // get the last one
  zavge /= zwt;
  vertex_list.push_back(zavge);    

  return;
}
*/

/*
// Used for non-INTT matching
void PHSiliconTpcTrackMatching::addSiliconClusters(  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>> &vertex_map)
{

  for(auto it = vertex_map.begin(); it != vertex_map.end(); ++it)
    {
      unsigned int tpcid = it->second.first;
      SvtxTrack *tpc_track = _track_map->get(tpcid);
      if(Verbosity() > 1) std::cout << "  tpcid " << tpcid << " original z " << tpc_track->get_z() << std::endl;
      
      // add the silicon cluster keys to the track
      unsigned int si_id = it->second.second;
      SvtxTrack *si_track = _track_map_silicon->get(si_id);
      if(Verbosity() > 1) std::cout << "  si track id " << si_id << std::endl;
      for (SvtxTrack::ConstClusterKeyIter si_iter = si_track->begin_cluster_keys();
	   si_iter != si_track->end_cluster_keys();
	   ++si_iter)
	{
	  TrkrDefs::cluskey si_cluster_key = *si_iter;
	  
	  if(Verbosity() > 1) 
	    cout << "   inserting si cluster key " << si_cluster_key << " into existing TPC track " << tpc_track->get_id() << endl;
	  
	  tpc_track->insert_cluster_key(si_cluster_key);
	  
	}
  
      // update the track position to the si one
      tpc_track->set_x(si_track->get_x());
      tpc_track->set_y(si_track->get_y());
      tpc_track->set_z(si_track->get_z());
      
      if(Verbosity() > 2)
	{
	  std::cout << " TPC seed track ID " << tpc_track->get_id() << " si track id " << si_track->get_id()
		    << " new nclus " << tpc_track->size_cluster_keys() << std::endl;
	}
      
      if(Verbosity() > 2)  
	tpc_track->identify();
    }
  
  return;
}	  
*/

/*
void PHSiliconTpcTrackMatching::cleanVertexMap(  
						   std::map<unsigned int, short int> &vertex_crossings_map,
						   std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map,
						   std::map<unsigned int, short int> &tpc_crossing_map )
{
  if(vertex_crossings_map.size() == 0) return;

  // Sanity check: is the initial rough crossing estimate consistent with the crossing number for this si vertex?

  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>> bad_map;
  for(auto [ivert, crossing] : vertex_crossings_map)
    {
      if(Verbosity() > 1) 
	std::cout << " CleanVertexMap:   ivert " << ivert << " crossing " << crossing << std::endl;		  
   
      // loop over all tracks associated with this vertex
      // remember that any tracks not associated with one of the vertices are unusable, they are ignored
      auto  ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  unsigned int tpcid = it->second.first;
	  unsigned int si_id = it->second.second;

	  // need the initial crossing estimate for comparison - stored in tpc_ crossing_map
	  auto sit = tpc_crossing_map.find(tpcid);	  
	  if(abs(sit->second - crossing) > 2)
	    {
	      // Fails the sanity check, mark for removal from the vertex map
	      if(Verbosity() > 1) 
		std::cout << "      Crossing mismatch: ivert " << ivert << " tpc id " << tpcid  << " vert crossing " << crossing << " track crossing " << sit->second << std::endl;

	      bad_map.insert(std::make_pair(ivert, std::make_pair(tpcid, si_id)));
	    }
	}
    }

  // If we delete an entry from vertex map, the TPC track is not associated with a bunch crossing, 
  //    -- the cluster z values are not changed, and the silicon clusters are not added. 
  // The track remains in the track map as a TPC-only track, by default assumed to be in bunch crossing zero
  // If multiple entries are deleted from vertex map, the TPC-only track copies are left in the track map, and also in _seed-track_map
  // The track cleaner will remove all but one of them, based on chisq from the Acts fitter.

  // remove bad entries from vertex_map so the wrong silicon is not associated
  for(auto [ivert, id_pair] : bad_map)
    {
      unsigned int tpc_id = id_pair.first;
      unsigned int si_id = id_pair.second;

      // Have to iterate over vertex_map and examine each pair to find the one matching bad_map
      // this logic works because we call the equal range on vertex_map for every id_pair
      // so we only delete one entry per equal range call
      auto ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  if(it->second.first == tpc_id && it->second.second == si_id)
	    {
	      if(Verbosity() > 1) 
		std::cout << "   erasing entry for ivert " << ivert << " tpc_id " << tpc_id << " si_id " << si_id << std::endl;
	      vertex_map.erase(it);
	      break;  // the iterator is no longer valid
	    }
	}
    } 
      
  return;
}
*/

/*
void PHSiliconTpcTrackMatching::getCrossingNumber(  
						  std::vector<double> &vertex_list,
						  std::multimap<unsigned int, std::pair<unsigned int, unsigned int>>  &vertex_map, 
						  std::map<unsigned int, short int> &vertex_crossings_map)
{
  if(vertex_list.size() == 0) return;

  for(unsigned int ivert = 0; ivert<vertex_list.size(); ++ivert)
    {     
      std::vector<double> crossing_vec;

      double vert_z = vertex_list[ivert];
      if(Verbosity() > 1) 
	std::cout << "Vertex " << ivert << " vertex z " << vert_z << std::endl;
      
      auto  ret = vertex_map.equal_range(ivert);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  unsigned int tpcid = it->second.first;
	  double tpc_z =  _track_map->get(tpcid)->get_z();
	  double z_mismatch = tpc_z - vert_z;

	  double crossings = getBunchCrossing(tpcid, z_mismatch);
	  crossing_vec.push_back(crossings);

	  if(Verbosity() > 1) 
	    std::cout << "   tpc_id " << tpcid << " tpc pT " << _track_map->get(tpcid)->get_pt() 
					<< " tpc eta " << _track_map->get(tpcid)->get_eta() << " si_id " << it->second.second 
					<< " tpc z " << tpc_z  << " crossings " << crossings << std::endl;
	}

      if(crossing_vec.size() == 0) continue;
      double crossing_median = getMedian(crossing_vec);
  
      double crossing_avge = 0.0;
      double crossing_wt = 0.0;
      for(auto cross : crossing_vec)
	{
	  if(fabs(cross-crossing_median) < 1.0)
	    {
	      crossing_avge += cross;
	      crossing_wt++;	      
	    }
	}
      short int crossing = SHRT_MAX;
      if(crossing_wt != 0)
	{
	  crossing_avge /= crossing_wt;      
	  double crossing_avge_rounded  =  round(crossing_avge);
	  crossing = (short int) crossing_avge_rounded;
	}
      vertex_crossings_map.insert(std::make_pair(ivert, crossing));

      if(Verbosity() > 1) 
	std::cout << "     crossing_median " << crossing_median << " crossing average = " << crossing_avge << " crossing_wt " << crossing_wt 
				    << " crossing integer " << crossing << std::endl;
    }
  return;
}
*/
/*
// Used for triggered crossing only
void PHSiliconTpcTrackMatching::addSiliconClusters( std::multimap<unsigned int, unsigned int> &tpc_matches )
{

  for(auto it = tpc_matches.begin(); it != tpc_matches.end(); ++it)
    {
      unsigned int tpcid = it->first;
      if(Verbosity() > 1) std::cout << "  tpcid " << tpcid << " original z " << _track_map->get(tpcid)->get_z() << std::endl;
      
      // add the silicon cluster keys to the track
      unsigned int si_id = it->second;
     
      auto svtxseed = std::make_unique<SvtxTrackSeed_v1>();
      svtxseed->set_silicon_seed_index(si_id);
      svtxseed->set_tpc_seed_index(tpcid);
      _svtx_seed_map->insert(svtxseed.get());

      if(Verbosity() > 1) std::cout << "  si track id " << si_id << std::endl;
     
    }

  return;
}	  
*/
/*
// used for INTT matching
void PHSiliconTpcTrackMatching::addSiliconClusters( std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches )
{

  for(auto it = crossing_matches.begin(); it != crossing_matches.end(); ++it)
    {
      unsigned int tpcid = it->second.first;
      if(Verbosity() > 1) std::cout << "  tpcid " << tpcid << " original z " << _track_map->get(tpcid)->get_z() << std::endl;
      
      // add the silicon cluster keys to the track
      unsigned int si_id = it->second.second;
      if(Verbosity() > 1) std::cout << "  si track id " << si_id << std::endl;
    
      auto seed = std::make_unique<SvtxTrackSeed_v1>();
      seed->set_silicon_seed_index(si_id);
      seed->set_tpc_seed_index(tpcid);
      _svtx_seed_map->insert(seed.get());
      
    }
  
  return;
}	  
*/

/*
void PHSiliconTpcTrackMatching::checkCrossingMatches( std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches,  
						      std::map<unsigned int, short int> &tpc_crossing_map )
{
  bool make_crossing_guess = false;

  float vdrift = 8.0e-3;
  //  float crossing_period = 106.0;

  std::multimap<short int, std::pair<unsigned int, unsigned int>> good_map;
  std::multimap<short int, std::pair<unsigned int, unsigned int>> bad_map;

  for(auto it = crossing_matches.begin(); it != crossing_matches.end(); ++it)
    {
      short int crossing = it->first;

      unsigned int tpcid = it->second.first;
      TrackSeed *tpc_track = _track_map->get(tpcid);
      
      unsigned int si_id = it->second.second;
      TrackSeed *si_track = _track_map_silicon->get(si_id);

      float tpc_eta = tpc_track->get_eta();
      float si_eta = si_track->get_eta();

      float z_si = si_track->get_z();
      float z_tpc = tpc_track->get_z();
      float z_mismatch = z_tpc-z_si;

      float mag_crossing_z_mismatch = fabs(crossing) * crossing_period * vdrift;

      // We do not know the sign  of the z mismatch for a given crossing unless we know the drift direction in the TPC, use mag
      // could instead look up any TPC cluster key in the track to get side
      // z-mismatch can occasionally be up to 2 crossings due to TPC extrapolation precision
      if( fabs( fabs(z_mismatch) - mag_crossing_z_mismatch ) < 3.0)
	{ 
	  if(Verbosity() > 1)	  
	    std::cout << "  Success:  crossing " << crossing << " tpc_eta " << tpc_eta << " si eta " << si_eta << " tpcid " << tpcid << " si id " << si_id 
		      << " tpc z " << z_tpc << " si z " << z_si << " z_mismatch " << z_mismatch << " mag_crossing_z_mismatch " << mag_crossing_z_mismatch << std::endl;
	}
      else
	{
	  if(Verbosity() > 1)
	    std::cout << "  FAILURE:  crossing " << crossing << " tpc_eta " << tpc_eta  << " si eta " << si_eta << " tpcid " << tpcid << " si id " << si_id 
		      << " tpc z " << z_tpc << " si z " << z_si << " z_mismatch " << z_mismatch << " mag_crossing_z_mismatch " << mag_crossing_z_mismatch << std::endl;

	  bad_map.insert(std::make_pair(crossing, std::make_pair(tpcid, si_id)));

	  if(make_crossing_guess)
	    {
	      // substitute a crossing estimate from the z_mismatch
	      short int crossing_guess = (short int) getBunchCrossing(tpcid, z_mismatch);
	      
	      if(Verbosity() > 1) std::cout << "         substitute crossing guess " << crossing_guess << " for crossing " << crossing << std::endl; 	      
	      good_map.insert(std::make_pair(crossing_guess, std::make_pair(tpcid, si_id)));
	    }
	}
    }

  // remove bad entries from crossing_matches
  for(auto [crossing, id_pair] : bad_map)
    {
      unsigned int tpcid = id_pair.first;
      unsigned int si_id = id_pair.second;
      
      // Have to iterate over crossing_matches and examine each pair to find the one matching bad_map
      // this logic works because we call the equal range on vertex_map for every id_pair
      // so we only delete one entry per equal range call
      auto ret = crossing_matches.equal_range(crossing);
      for(auto it = ret.first; it != ret.second; ++it)
	{
	  if(it->second.first == tpcid && it->second.second == si_id)
	    {
	      if(Verbosity() > 0) std::cout <<  "   checkCrossingMatches: erasing tpc_crossing_map entry for crossing " << crossing << " tpcid " << tpcid  << std::endl;
	      tpc_crossing_map.erase(tpcid);

	      if(Verbosity() > 1) 
		std::cout << "                        erasing crossing_matches entry for crossing " << crossing << " tpcid " << tpcid << " si_id " << si_id << std::endl;
	      crossing_matches.erase(it);
	      break;  // the iterator is no longer valid
	    }
	}
    }	  

  if(make_crossing_guess)
    {
      // replace them with crossing guess
      for(auto [crossing_guess, id_pair] : good_map)
	{
	  unsigned int tpcid = id_pair.first;
	  unsigned int si_id = id_pair.second;
	  
	  if(Verbosity() > 1) 
	    std::cout << "  checkCrossingMatches:  adding crossing_matches and tpc_crossing_map entry for crossing_guess " << crossing_guess << " tpcid " << tpcid 
		      << " si_id " << si_id << std::endl;
	  crossing_matches.insert(std::make_pair(crossing_guess,std::make_pair(tpcid, si_id)));
	  tpc_crossing_map.insert(std::make_pair(tpcid, crossing_guess));
	}
    }
    
  return;
}	  
*/

/*
// uses INTT time to get bunch crossing
void PHSiliconTpcTrackMatching::getMatchCrossingIntt(  
						std::multimap<unsigned int, unsigned int> &tpc_matches,
						std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches,
						std::map<unsigned int, short int> &tpc_crossing_map )
{
  if(Verbosity() > 0) std::cout << " tpc_matches size " << tpc_matches.size() << std::endl;

  for(auto [tpcid, si_id] : tpc_matches)
    {
      TrackSeed *si_track =_track_map_silicon->get(si_id);

      if(Verbosity() > 0) 
	std::cout << "TPC track " << tpcid << " si track " << si_id << std::endl;

      // If the Si track contains an INTT hit, use it to get the bunch crossing offset

      std::vector<short int> intt_crossings = getInttCrossings(si_track);      

      bool keep_it = true;
      short int crossing_keep = 0;
      if(intt_crossings.size() == 0) 
	{
	  if(Verbosity() > 0)  std::cout << " Silicon track " << si_id << " has no INTT clusters, skip this combination " << std::endl;
	  continue ;
	}
      else
	{
	  crossing_keep = intt_crossings[0];
	  for(unsigned int ic=1; ic<intt_crossings.size(); ++ic)
	    {	  
	      if(intt_crossings[ic] != crossing_keep)
		{
		  if(Verbosity() > -1)
		    {
		      std::cout << " Warning: INTT crossings not all the same for tpc track " << tpcid << " silicon track " << si_id << " crossing_keep " << crossing_keep << " new crossing " << intt_crossings[ic] << " keep the first one in the list" << std::endl;  
		    }
		}
	    }
	}
      
      if(keep_it)
	{            
	  crossing_matches.insert(std::make_pair(crossing_keep,std::make_pair(tpcid, si_id)));
	  tpc_crossing_map.insert(std::make_pair(tpcid, crossing_keep));
	  
	  if(Verbosity() > 0) 
	    std::cout << "                    tpc track " << tpcid << " si Track " << si_id << " final crossing " << crossing_keep  << std::endl;           
	}
    }
  return;
}
*/

/*
void PHSiliconTpcTrackMatching::addTrackBunchCrossing( std::map<unsigned int, short int> &tpc_crossing_map)
{
  if(tpc_crossing_map.size() == 0) return;
  
  for(auto [tpcid, crossing] : tpc_crossing_map)
    {
     if(Verbosity() > 1) 
      std::cout << PHWHERE << "  Add bunch crossing to track " << tpcid << " crossing " << crossing << std::endl;		  

       TrackSeed *tpc_track = _track_map->get(tpcid);
       tpc_track->set_crossing(crossing);
    }
  return;
}
*/
/*
double PHSiliconTpcTrackMatching::getMedian(std::vector<double> &v)
{
  if(v.size() == 0) return NAN;

  double median = 0.0;

  if( (v.size() % 2) == 0)
    {
      // even number of entries
      // we want the average of the middle two numbers, v.size()/2 and v.size()/2-1
      auto m1 = v.begin() + v.size()/2;
      std::nth_element(v.begin(), m1, v.end());
      double median1 =  v[v.size()/2]; 

      auto m2 = v.begin() + v.size()/2 - 1;
      std::nth_element(v.begin(), m2, v.end());
      double median2 =  v[v.size()/2 - 1]; 

      median = (median1 + median2) / 2.0; 
      if(Verbosity() > 2) std::cout << "The vector size is " << v.size() 
				    << " element m is " << v.size() / 2  << " = " << v[v.size()/2] 
				    << " element m-1 is " << v.size() / 2 -1 << " = " << v[v.size()/2-1] 
				    <<  std::endl;
    } 
  else
    {
      // odd number of entries
      auto m = v.begin() + v.size()/2;
      std::nth_element(v.begin(), m, v.end());
      median =  v[v.size()/2];
      if(Verbosity() > 2) std::cout << "The vector size is " << v.size() << " element m is " << v.size() / 2 << " = " << v[v.size()/2] <<  std::endl;
    }

    return median ;
}
*/
/*
double PHSiliconTpcTrackMatching::getBunchCrossing(unsigned int trid, double z_mismatch )
{
  double vdrift = 8.00;  // cm /microsecond
  //double z_bunch_separation = 0.106 * vdrift;  // 106 ns bunch crossing interval, as in pileup generator
  double z_bunch_separation = (crossing_period/1000.0) * vdrift;  // 106 ns bunch crossing interval, as in pileup generator

  // The sign of z_mismatch will depend on which side of the TPC the tracklet is in
  TrackSeed *track = _track_map->get(trid);

  double crossings = z_mismatch / z_bunch_separation;

  // Check the TPC side for the first cluster in the track
  unsigned int side = 10;
  for (TrackSeed::ConstClusterKeyIter iter = track->begin_cluster_keys();
       iter != track->end_cluster_keys();
       ++iter)
    {
      TrkrDefs::cluskey cluster_key = *iter;
      unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
      if(trkrid == TrkrDefs::tpcId)
	{
	  side = TpcDefs::getSide(cluster_key);
	  break;   // we only need the first one  
	}
    }

  if(side == 10) return SHRT_MAX;

  // if side = 1 (north, +ve z side), a positive t0 will make the cluster late relative to true z, so it will look like z is less positive
  // so a negative z mismatch for side 1 means a positive t0, and positive crossing, so reverse the sign for side 1
  if(side == 1)
    crossings *= -1.0;
  
  if(Verbosity() > 1) 
    std::cout << "             trackid " << trid << " side " << side << " z_mismatch " << z_mismatch << " crossings " << crossings << std::endl;

  return crossings;
}
*/
/*
void PHSiliconTpcTrackMatching::tagInTimeTracks(  
						std::multimap<unsigned int, unsigned int> &tpc_matches,
						std::multimap<int, std::pair<unsigned int, unsigned int>> &crossing_matches,
						std::map<unsigned int, int> &tpc_crossing_map )
{

  for(auto [tpcid, si_id] : tpc_matches)
    {
      TrackSeed *tpc_track = _track_map->get(tpcid);
      double tpc_z = tpc_track->get_z();
      double tpc_pt = tpc_track->get_pt();

      TrackSeed *si_track =_track_map_silicon->get(si_id);
      double si_z = si_track->get_z();

      // this factor will increase the window size at low pT
      // otherwise the matching efficiency drops off at low pT
      // it would be better if this was a smooth function
      double mag = 1.0;
      if(tpc_pt < 6.0) mag = 2;
      if(tpc_pt < 3.0)  mag = 4.0;

      // Check for a match in z
      if(fabs(tpc_z - si_z) < _z_search_win * mag)
	{
	  //track from  triggered event
	  int crossing = 0;
	  crossing_matches.insert(std::make_pair(crossing,std::make_pair(tpcid, si_id)));
	  tpc_crossing_map.insert(std::make_pair(tpcid, crossing));
	  if(Verbosity() > 1)
	    std::cout << " triggered: tpc_trackid " << tpcid
		      << " eta " << tpc_track->get_eta()
		      << " pT " << tpc_track->get_pt() 
		      << " si_trackid " << si_id 
		      << " z_tpc " << tpc_track->get_z() 
		      << " z_si " << si_track->get_z() 
		      << " crossing " << crossing 
		      << std::endl;				  
	}
    }
  return;
}

void PHSiliconTpcTrackMatching::tagMatchCrossing(  
						std::multimap<unsigned int, unsigned int> &tpc_matches,
						std::multimap<short int, std::pair<unsigned int, unsigned int>> &crossing_matches,
						std::map<unsigned int, short int> &tpc_crossing_map )
{

  for(auto [tpcid, si_id] : tpc_matches)
    {
      TrackSeed *tpc_track = _track_map->get(tpcid);
      double tpc_z = tpc_track->get_z();

      TrackSeed *si_track =_track_map_silicon->get(si_id);
      double si_z = si_track->get_z();

      // this is an initial estimate of the bunch crossing based on the z-mismatch for this track
      short int crossing = (short int) getBunchCrossing(tpcid, tpc_z - si_z);
      crossing_matches.insert(std::make_pair(crossing,std::make_pair(tpcid, si_id)));
      tpc_crossing_map.insert(std::make_pair(tpcid, crossing));
      if(Verbosity() > 1)
	std::cout << " pileup: tpc_trackid " << tpcid
		  << " eta " << tpc_track->get_eta()
		  << " pT " << tpc_track->get_pt() 
		  << " si_trackid " << si_id
		  << " z_tpc " << tpc_track->get_z() 
		  << " z_si " << si_track->get_z() 
		  << " crossing " << crossing 
		  << std::endl;		
    }

return;
}

// This is for non-pp mode, i.e. straight geometric matching including a z cut
void PHSiliconTpcTrackMatching::addTrackBunchCrossing(std::multimap<unsigned int, unsigned int> &tpc_matches)	  
{
  if(tpc_matches.size() == 0) return;
  
  for(auto [tpcid, si_id] : tpc_matches)
    {
      short int crossing = 0;

       TrackSeed *tpc_track = _track_map->get(tpcid);
       TrackSeed *si_track = _track_map_silicon->get(si_id);
       if(!tpc_track)
	 {
	   std::cout << PHWHERE << "Did not find track " << tpcid << std::endl;
	   continue;
	 }

       if(Verbosity() > 1) 
	 std::cout << PHWHERE << "  Add bunch crossing to track " << tpcid << " crossing " << crossing << std::endl;		  
       
       si_track->set_crossing(crossing);
       tpc_track->set_crossing(crossing);	 
    }

  return;
}
*/
