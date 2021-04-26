#include "PHSiliconTpcTrackMatching.h"

#include "AssocInfoContainer.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

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
 PHTrackPropagating(name)
 , _track_map_name_silicon("SvtxSiliconTrackMap")
{
  //cout << "PHSiliconTpcTrackMatching::PHSiliconTpcTrackMatching(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHSiliconTpcTrackMatching::~PHSiliconTpcTrackMatching()
{

}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::Setup(PHCompositeNode *topNode)
{
  // put these in the output file
  cout << PHWHERE "_is_ca_seeder " << _is_ca_seeder << " Search windows: phi " << _phi_search_win << " eta " 
       << _eta_search_win << endl;

  // corrects the PHTpcTracker phi bias
  fdphi = new TF1("f1", "[0] + [1]/x^[2]");
  fdphi->SetParameter(0, _par0);
  fdphi->SetParameter(1, _par1);
  fdphi->SetParameter(2, _par2);

  // corrects the space charge distortion phi bias
  if(!_is_ca_seeder)
    {
      // PHTpcTracker correction is opposite in sign
      // and different in magnitude - why?
      _parsc0 *= -1.0 * 0.7;
      _parsc1 *= -1.0 * 0.7;
    }
  fscdphi = new TF1("f2","[0] + [1]*x^2");
  fscdphi->SetParameter(0, _parsc0 * _collision_rate / _reference_collision_rate);
  fscdphi->SetParameter(1, _parsc1 * _collision_rate / _reference_collision_rate);

  int ret = PHTrackPropagating::Setup(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::Process()
{
  // _track_map contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // We will add the silicon clusters to the TPC tracks already on the node tree
  // We will have to expand the number of tracks whenever we find multiple matches to the silicon

  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size() << " Silicon track map size " << _track_map_silicon->size() << endl;

  if(Verbosity() > 2)
    {
      // list silicon tracks
      for (auto phtrk_iter_si = _track_map_silicon->begin();
	   phtrk_iter_si != _track_map_silicon->end(); 
	   ++phtrk_iter_si)
	{
	  _tracklet_si = phtrk_iter_si->second;	  
	  
	  double si_phi = atan2(_tracklet_si->get_py(), _tracklet_si->get_px());
	  double si_eta = _tracklet_si->get_eta();
	  
	  cout << " Si track " << _tracklet_si->get_id()  << " si_phi " << si_phi  << " si_eta " << si_eta << endl;
	}  
    }
  
  
  // We remember the original size of the TPC track map here
  const unsigned int original_track_map_lastkey = _track_map->end()->first;
  
  // loop over the original TPC tracks
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first >= original_track_map_lastkey)  break;
      
      _tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet_tpc->get_phi()
	    << endl;
	}

      double tpc_phi = atan2(_tracklet_tpc->get_py(), _tracklet_tpc->get_px());
      double tpc_eta = _tracklet_tpc->get_eta();
      double tpc_pt = sqrt( pow(_tracklet_tpc->get_px(),2) + pow(_tracklet_tpc->get_py(),2) );

      // phi correction for PHTpcTracker tracklets is charge dependent
      double sign_phi_correction = _tracklet_tpc->get_charge();

      /// Correct the correction for the field direction
      /// Kludge to get the phi matching correct based on the field
      /// direction
      if(_field.find(".root") != std::string::npos)
	sign_phi_correction *= -1;
      if(_fieldDir > 0)
	sign_phi_correction *= -1;

      // hard code this here for now
      // this factor will increase the window size at low pT
      // otherwise the matching efficiency drops off at low pT
      // it would be better if this was a smooth function
      double mag = 1.0;
      if(tpc_pt < 6.0) mag = 2;
      if(tpc_pt < 3.0)  mag = 4.0;

      if(Verbosity() > 3)
	{
	  cout << "Original TPC tracklet:" << endl;
	  _tracklet_tpc->identify();
	}

      // correct the TPC tracklet phi for the space charge offset, if this is the calib pass
      // this is done just to let us tighten up the matching window

      if(_sc_calib_flag)
	{
	  tpc_phi -= fscdphi->Eval(tpc_eta);
	}
      // the distortion correction can push tpc_phi outside +/- M_PI
      if(tpc_phi < - M_PI) tpc_phi += 2.0*M_PI;
      if(tpc_phi > M_PI) tpc_phi -= 2.0*M_PI;


      // Now search the silicon track list for a match in eta and phi
      // NOTE: we will take the combined track vertex from the vertex associated with the silicon stub, once the match is made

      std::set<unsigned int> si_matches;
      for (auto phtrk_iter_si = _track_map_silicon->begin();
	   phtrk_iter_si != _track_map_silicon->end(); 
	   ++phtrk_iter_si)
	{
	  _tracklet_si = phtrk_iter_si->second;	  

	  double si_phi = atan2(_tracklet_si->get_py(), _tracklet_si->get_px());
	  double si_eta = _tracklet_si->get_eta();
	  //double si_pt = sqrt(pow(_tracklet_si->get_px(), 2) + pow(_tracklet_si->get_py(), 2) );

	  if(Verbosity() >= 2)
	    {
	      cout << " testing for a match for TPC track " << _tracklet_tpc->get_id() << " with Si track " << _tracklet_si->get_id() << endl;	  
	      cout << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " <<   tpc_phi-si_phi << " phi search " << _phi_search_win  << " tpc_eta " << tpc_eta 
		   << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << " eta search " << _eta_search_win << endl;
	    }

	  bool eta_match = false;
	  bool phi_match = false;
	  if(  fabs(tpc_eta - si_eta) < _eta_search_win * mag) eta_match = true;

	  // PHTpcTracker has a bias in the tracklet phi that depends on charge sign, PHCASeeding does not
	  if(_is_ca_seeder)
	    {
	      if(  fabs(tpc_phi - si_phi)  < _phi_search_win * mag) phi_match = true;
	    }
	  else
	    {
	      // PHTpcTracker
	      //double si_pt = sqrt( pow(_tracklet_si->get_px(),2) + pow(_tracklet_si->get_py(),2) );
	      double phi_search_win_lo = fdphi->Eval(tpc_pt) * sign_phi_correction -  _phi_search_win * mag;
	      double phi_search_win_hi = fdphi->Eval(tpc_pt) * sign_phi_correction +  _phi_search_win * mag;

	      if(Verbosity() > 10) 
		cout << " phi_search_win_lo " << phi_search_win_lo << " phi_search_win_hi " << phi_search_win_hi << endl;

	      if(  (tpc_phi - si_phi) > phi_search_win_lo && (tpc_phi - si_phi) < phi_search_win_hi) phi_match = true;	      
	    }

	  /*
	  // temporary for debugging!
	  if(_test_windows)
	    cout << " Try_silicon:  pt " << tpc_pt << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi-si_phi  
		 << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << endl;
	  */	  
	  
	  if(eta_match && phi_match)
	    {
	      // got a match, add to the list
	      if(Verbosity() >= 1)  
		{
		  cout << " found a match for TPC track " << _tracklet_tpc->get_id() << " with Si track " << _tracklet_si->get_id() << endl;
		  cout << "          tpc_phi " << tpc_phi << " si_phi " <<  si_phi << " phi_match " << phi_match 
		       << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " eta_match " << eta_match << endl;
		}

	      // temporary!
	      if(_test_windows)
		cout << " Try_silicon:  pt " << tpc_pt << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi-si_phi  
		     << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << endl;

	      si_matches.insert(_tracklet_si->get_id());
	    }
	  else
	    {
	      if(Verbosity() >= 10) 
		{ 
		  cout << " no match for TPC track " << _tracklet_tpc->get_id() << " with Si track " << _tracklet_si->get_id() << endl;
		  cout << "          tpc_phi " << tpc_phi << " si_phi " <<  si_phi << " phi_match " << phi_match 
		       << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " eta_match " << eta_match << endl;
		}	  
	    }
	}
      // we did not get a match, sound the alarm
      if(si_matches.size() == 0)
	{
	  if(Verbosity() >= 1)
	    {
	      cout << PHWHERE << " Did NOT find a match for TPC track " << _tracklet_tpc->get_id()  << "  tpc_phi " << tpc_phi << " tpc_eta " << tpc_eta  << endl;
	    }

	  // we did not get a silicon track match, so set the track vertex arbitrarily to vertex 0 if one does not exist already
	  unsigned int vertexId = _tracklet_tpc->get_vertex_id();
	  if(vertexId == UINT_MAX)  vertexId = 0;
	  _tracklet_tpc->set_vertex_id(vertexId);
	  
	  // set the track position to the vertex position
	  const SvtxVertex *svtxVertex = _vertex_map->get(vertexId);      
	  if(svtxVertex)
	    {
	      _tracklet_tpc->set_x(svtxVertex->get_x());
	      _tracklet_tpc->set_y(svtxVertex->get_y());
	      _tracklet_tpc->set_z(svtxVertex->get_z());
	    }
	  else
	    {
	      std::cout << PHWHERE << "Failed to find vertex object for vertex ID " << vertexId 
			<< " associated with TPC tracklet " << _tracklet_tpc->get_id() 
			<< " si tracklet " << _tracklet_si->get_id()
			<< " set vertex to (0,0,0)"
			<< std::endl; 
	      std::cout << " ---------- Silicon track --------------" << std::endl;
	      _tracklet_si->identify();
	      std::cout << " ---------- TPC track -----------------" << std::endl;
	      _tracklet_tpc->identify();

	      _tracklet_tpc->set_x(0.0);
	      _tracklet_tpc->set_y(0.0);
	      _tracklet_tpc->set_z(0.0);
	    }
	}

      // Add the silicon clusters to the track
      unsigned int isi = 0;
      for(auto si_it = si_matches.begin(); si_it != si_matches.end(); ++si_it)
	{
	  // get the si tracklet for this id
	  _tracklet_si = _track_map_silicon->get(*si_it);

	  if(Verbosity() > 0)
	    cout << "   isi = " << isi << " si tracklet " << _tracklet_si->get_id() << " was matched to TPC tracklet " << _tracklet_tpc->get_id() << endl;

	  // get the silicon clusters
	  std::set<TrkrDefs::cluskey> si_clusters;
	  for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_si->begin_cluster_keys();
	       key_iter != _tracklet_si->end_cluster_keys();
	       ++key_iter)
	    {
	      TrkrDefs::cluskey cluster_key = *key_iter;
	      if(Verbosity() >=1) cout << "   inserting cluster key " << cluster_key << " into list " << " for si tracklet " << _tracklet_si->get_id() << endl;
	      si_clusters.insert(cluster_key);
	    }

	  if(isi == 0)
	    {
	      // update the original track on the node tree

	      // the track takes its vertex from the si stub
	      unsigned int vertexId = _tracklet_si->get_vertex_id();
	      _tracklet_tpc->set_vertex_id(vertexId);

	      // set the track position to the vertex position
	      const SvtxVertex *svtxVertex = _vertex_map->get(vertexId);      
	      if(svtxVertex)
		{
		  _tracklet_tpc->set_x(svtxVertex->get_x());
		  _tracklet_tpc->set_y(svtxVertex->get_y());
		  _tracklet_tpc->set_z(svtxVertex->get_z());
		}
	      else
		{
		  std::cout << PHWHERE << "Failed to find vertex object for vertex ID " << vertexId 
			    << " associated with TPC tracklet " << _tracklet_tpc->get_id() 
			    << " si tracklet " << _tracklet_si->get_id()
			    << " set vertex to (0,0,0)"
			    << std::endl; 
		  std::cout << " ---------- Silicon track --------------" << std::endl;
		  _tracklet_si->identify();
		  std::cout << " ---------- TPC track -----------------" << std::endl;
		  _tracklet_tpc->identify();

		  _tracklet_tpc->set_x(0.0);
		  _tracklet_tpc->set_y(0.0);
		  _tracklet_tpc->set_z(0.0);
		}
	      for(auto clus_iter=si_clusters.begin(); clus_iter != si_clusters.end(); ++clus_iter)
		{
		  if(Verbosity() >= 1) cout << "   inserting si cluster key " << *clus_iter << " into exisiting TPC track " << _tracklet_tpc->get_id() << endl;
		  _tracklet_tpc->insert_cluster_key(*clus_iter);
		  _assoc_container->SetClusterTrackAssoc(*clus_iter, _tracklet_tpc->get_id());
		}

	      if(Verbosity() > 3)
		_tracklet_tpc->identify();
	    }
	  else
	    {
	      // more than one si stub matches
	      // make a copy of the TPC track, update it and add it to the end of the node tree 
	      #if __cplusplus < 201402L
	      auto newTrack = boost::make_unique<SvtxTrack_v2>();
	      #else
	      auto newTrack = std::make_unique<SvtxTrack_v2>();
	      #endif
	      const unsigned int lastTrackKey = _track_map->end()->first; 
	      if(Verbosity() >= 1) cout << "Extra match, add a new track to node tree with key " <<  lastTrackKey << endl;
	      
	      newTrack->set_id(lastTrackKey);

	      unsigned int vertexId = _tracklet_si->get_vertex_id();
	      newTrack->set_vertex_id(vertexId);

	      // set the track position to the vertex position
	      const SvtxVertex *svtxVertex = _vertex_map->get(vertexId);      
	      if(svtxVertex)
		{
		  newTrack->set_x(svtxVertex->get_x());
		  newTrack->set_y(svtxVertex->get_y());
		  newTrack->set_z(svtxVertex->get_z());
		}
	      else
		{
		  std::cout << PHWHERE << "Failed to find vertex object for vertex ID " << vertexId << " associated with TPC tracklet " << _tracklet_tpc->get_id() << std::endl; 
		  std::cout << " ---------- Silicon track --------------" << std::endl;
		  _tracklet_si->identify();
		  std::cout << " ---------- TPC track -----------------" << std::endl;
		  _tracklet_tpc->identify();
		  newTrack->set_x(0.0);
		  newTrack->set_y(0.0);
		  newTrack->set_z(0.0);
		}
		      
	      newTrack->set_charge(_tracklet_tpc->get_charge());
	      newTrack->set_px(_tracklet_tpc->get_px());
	      newTrack->set_py(_tracklet_tpc->get_py());
	      newTrack->set_pz(_tracklet_tpc->get_pz());
	      newTrack->set_x(_tracklet_tpc->get_x());
	      newTrack->set_y(_tracklet_tpc->get_y());
	      newTrack->set_z(_tracklet_tpc->get_z());
	      for(int i = 0; i < 6; ++i)
		{
		  for(int j = 0; j < 6; ++j)
		    {
		      newTrack->set_error(i,j, _tracklet_tpc->get_error(i,j));
		    }
		}
	      
	      // loop over associated clusters to get hits for TPC only, add to new track copy
	      for (SvtxTrack::ConstClusterKeyIter iter = _tracklet_tpc->begin_cluster_keys();
		   iter != _tracklet_tpc->end_cluster_keys();
		   ++iter)
		{
		  TrkrDefs::cluskey cluster_key = *iter;
		  unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
		  if(trkrid == TrkrDefs::tpcId)
		    {
		      newTrack->insert_cluster_key(cluster_key);
		      _assoc_container->SetClusterTrackAssoc(cluster_key, _tracklet_tpc->get_id());
		    }
		}
		  
	      // now add the new si clusters
	      for(auto clus_iter=si_clusters.begin(); clus_iter != si_clusters.end(); ++clus_iter)
		{
		  if(Verbosity() >= 1) cout << "   inserting si cluster key " << *clus_iter << " into new track " << newTrack->get_id() << endl;
		  newTrack->insert_cluster_key(*clus_iter);
		  _assoc_container->SetClusterTrackAssoc(*clus_iter, newTrack->get_id());
		}

	      if(Verbosity() > 3)
		{
		  cout << "  -- inserting new track with id " << newTrack->get_id() << " into trackmap " << endl;
		  newTrack->identify();
		}

	      _track_map->insert(newTrack.get());
	 
	    }
	  
	  isi++;
	}
    }

  if(Verbosity() > 0)  
    cout << " Final track map size " << _track_map->size() << endl;
  
  if (Verbosity() >= 1)
    cout << "PHSiliconTpcTrackMatching::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconTpcTrackMatching::End()
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int  PHSiliconTpcTrackMatching::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get additional objects off the Node Tree
  //---------------------------------

  _track_map_silicon = findNode::getClass<SvtxTrackMap>(topNode,  "SvtxSiliconTrackMap");
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxSiliconTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

