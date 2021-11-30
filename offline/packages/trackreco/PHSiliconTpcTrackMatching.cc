#include "PHSiliconTpcTrackMatching.h"

#include "AssocInfoContainer.h"

/// Tracking includes
#include <trackbase/TrkrDefs.h>                // for cluskey, getTrkrId, tpcId
#include <trackbase/TpcSeedTrackMapv1.h>     
#include <trackbase/TrkrClusterContainerv3.h>   
#include <trackbase/TrkrClusterv3.h>   
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitDefs.h>  // for keytype

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>


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
  SubsysReco(name)
 , _track_map_name_silicon("SvtxSiliconTrackMap")
{
  //cout << "PHSiliconTpcTrackMatching::PHSiliconTpcTrackMatching(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHSiliconTpcTrackMatching::~PHSiliconTpcTrackMatching()
{

}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::InitRun(PHCompositeNode *topNode)
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
 
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return ret;
}

//____________________________________________________________________________..
int PHSiliconTpcTrackMatching::process_event(PHCompositeNode*)
{
  // _track_map contains the TPC seed track stubs
  // _track_map_silicon contains the silicon seed track stubs
  // We will add the silicon clusters to the TPC tracks already on the node tree
  // We will have to expand the number of tracks whenever we find multiple matches to the silicon

  _seed_track_map->Reset();

  if(Verbosity() > 0)
    cout << PHWHERE << " TPC track map size " << _track_map->size() << " Silicon track map size " << _track_map_silicon->size() << endl;

  if(_track_map->size() == 0)
    return Fun4AllReturnCodes::EVENT_OK;
    
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
  const unsigned int original_track_map_lastkey = _track_map->empty() ? 0:std::prev(_track_map->end())->first;
  if(Verbosity() > 1)
    std::cout << "Original track map has lastkey " << original_track_map_lastkey << std::endl;

  // loop over the original TPC tracks
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first > original_track_map_lastkey)  break;
      
      _tracklet_tpc = phtrk_iter->second;
      
      if (Verbosity() > 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet_tpc-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet_tpc->get_phi()
	    << endl;
	}

      // This will always end up as a final track, no matter what - add it to the seed-track-map
      if(Verbosity() > 1) 
	std::cout << " TPC seed ID " << _tracklet_tpc->get_id() << " original nclus " << _tracklet_tpc->size_cluster_keys() << std::endl;
      _seed_track_map->addAssoc(_tracklet_tpc->get_id(), _tracklet_tpc->get_id());

      double tpc_phi = atan2(_tracklet_tpc->get_py(), _tracklet_tpc->get_px());
      double tpc_eta = _tracklet_tpc->get_eta();
      double tpc_pt = sqrt( pow(_tracklet_tpc->get_px(),2) + pow(_tracklet_tpc->get_py(),2) );

      // phi correction for PHTpcTracker tracklets is charge dependent
      double sign_phi_correction = _tracklet_tpc->get_charge();

      /// Correct the correction for the field direction
      /// Kludge to get the phi matching correct based on the field
      /// direction
      if(_field.find("2d") != std::string::npos)
	{
	  sign_phi_correction *= -1;
	  if(_fieldDir > 0)
	    sign_phi_correction *= -1;
	}
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

      double tpc_x = _tracklet_tpc->get_x();
      double tpc_y = _tracklet_tpc->get_y();
      double tpc_z = _tracklet_tpc->get_z();

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
	  double si_x = _tracklet_si->get_x();
	  double si_y = _tracklet_si->get_y();
	  double si_z = _tracklet_si->get_z();

	  if(Verbosity() >= 2)
	    {
	      cout << " testing for a match for TPC track " << _tracklet_tpc->get_id() << " with Si track " << _tracklet_si->get_id() << endl;	  
	      cout << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " <<   tpc_phi-si_phi << " phi search " << _phi_search_win  << " tpc_eta " << tpc_eta 
		   << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << " eta search " << _eta_search_win << endl;
	      std::cout << "      tpc x " << tpc_x << " si x " << si_x << " tpc y " << tpc_y << " si y " << si_y << " tpc_z " << tpc_z  << " si z " << si_z << std::endl;
	    }

	  bool eta_match = false;
	  bool phi_match = false;
	  bool position_match = false;
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

	  if(
	     fabs(tpc_x - si_x) < _x_search_win * mag &&
	     fabs(tpc_y - si_y) < _y_search_win * mag &&
	     fabs(tpc_z - si_z) < _z_search_win * mag
	     )
	    position_match = true;
	    
	  /*
	  // temporary for debugging!
	  if(_test_windows)
	    cout << " Try_silicon:  pt " << tpc_pt << " tpc_phi " << tpc_phi << " si_phi " << si_phi << " dphi " << tpc_phi-si_phi  
		 << " tpc_eta " << tpc_eta << " si_eta " << si_eta << " deta " << tpc_eta-si_eta << endl;
	  */	  
	  
	  if(eta_match && phi_match && position_match)
	    {
	      // got a match, add to the list
	      if(Verbosity() > 1)  
		{
		  cout << " found a match for TPC track " << _tracklet_tpc->get_id() << " with Si track " << _tracklet_si->get_id() << endl;
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
	  if(Verbosity() > 1)
	    {
	      cout << PHWHERE << " Did NOT find a match for TPC track " << _tracklet_tpc->get_id()  << "  tpc_phi " << tpc_phi << " tpc_eta " << tpc_eta  << endl;
	    }

	  // The TPC track seed vertex association is already done
	  // The TPC track seed position has already been set in the vertex associator to the PCA from the helical fit to clusters
	  // we do not need to change them here
	}
	  
      // Add the silicon clusters to the track if there is one or more matches
      unsigned int isi = 0;
      for(auto si_it = si_matches.begin(); si_it != si_matches.end(); ++si_it)
	{
	  // get the si tracklet for this id
	  _tracklet_si = _track_map_silicon->get(*si_it);
	  
	  if(Verbosity() > 1)
	    cout << "   isi = " << isi << " si tracklet " << _tracklet_si->get_id() << " was matched to TPC tracklet " << _tracklet_tpc->get_id() << endl;
	  
	  // get the silicon clusters
	  std::set<TrkrDefs::cluskey> si_clusters;
	  for (SvtxTrack::ConstClusterKeyIter key_iter = _tracklet_si->begin_cluster_keys();
	       key_iter != _tracklet_si->end_cluster_keys();
	       ++key_iter)
	    {
	      TrkrDefs::cluskey cluster_key = *key_iter;
	      if(Verbosity() >1) cout << "   inserting cluster key " << cluster_key << " into list " << " for si tracklet " << _tracklet_si->get_id() << endl;
	      si_clusters.insert(cluster_key);
	    }
	  
	  if(isi == 0)
	    {
	      // update the original track on the node tree
	      _tracklet_tpc->set_x(_tracklet_si->get_x());
	      _tracklet_tpc->set_y(_tracklet_si->get_y());
	      _tracklet_tpc->set_z(_tracklet_si->get_z());
	      
	      for(auto clus_iter=si_clusters.begin(); clus_iter != si_clusters.end(); ++clus_iter)
		{
		  if(Verbosity() > 1) 
		    cout << "   inserting si cluster key " << *clus_iter << " into existing TPC track " << _tracklet_tpc->get_id() << endl;
		  
		  _tracklet_tpc->insert_cluster_key(*clus_iter);
		  _assoc_container->SetClusterTrackAssoc(*clus_iter, _tracklet_tpc->get_id());
		}
	      
	      if(Verbosity() > 1)
		std::cout << " TPC seed track ID " << _tracklet_tpc->get_id() 
			  << " new nclus " << _tracklet_tpc->size_cluster_keys() << std::endl;
	      
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
	      const unsigned int lastTrackKey =  _track_map->empty() ? 0:std::prev(_track_map->end())->first; 
	      if(Verbosity() > 1) cout << "Extra match, add a new track to node tree with key " <<  lastTrackKey + 1 << endl;
	      
	      newTrack->set_id(lastTrackKey+1);
	    
	      newTrack->set_x(_tracklet_si->get_x());
	      newTrack->set_y(_tracklet_si->get_y());
	      newTrack->set_z(_tracklet_si->get_z());
	      
	      newTrack->set_charge(_tracklet_tpc->get_charge());
	      newTrack->set_px(_tracklet_tpc->get_px());
	      newTrack->set_py(_tracklet_tpc->get_py());
	      newTrack->set_pz(_tracklet_tpc->get_pz());
	      
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
		  if(Verbosity() > 1) 
		    cout << "   inserting si cluster key " << *clus_iter << " into new track " << newTrack->get_id() << endl;
		  newTrack->insert_cluster_key(*clus_iter);
		  _assoc_container->SetClusterTrackAssoc(*clus_iter, newTrack->get_id());
		}
	      
	      if(Verbosity() > 3)
		{
		  cout << "  -- inserting new track with id " << newTrack->get_id() << " into trackmap " << endl;
		  newTrack->identify();
		}
	      
	      _track_map->insert(newTrack.get());
	      
	      if(Verbosity() > 1)
		std::cout << " TPC seed track ID " << _tracklet_tpc->get_id() 
			  << " new track ID " << newTrack->get_id() << " new nclus " << newTrack->size_cluster_keys() << std::endl;
	      _seed_track_map->addAssoc(_tracklet_tpc->get_id(), newTrack->get_id() ) ;	 
	    }
	  
	  isi++;
	}
    }

  // loop over all tracks and copy the silicon clusters to the corrected cluster map
  if(_corrected_cluster_map)
    copySiliconClustersToCorrectedMap();

  
  if(Verbosity() > 0)  
    cout << " Final track map size " << _track_map->size() << " seed-track map size " << _seed_track_map->size() << endl;
  
  if (Verbosity() >= 1)
    cout << "PHSiliconTpcTrackMatching::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
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

   _assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
  if (!_assoc_container)
  {
    cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map_silicon = findNode::getClass<SvtxTrackMap>(topNode, _silicon_track_map_name);
  if (!_track_map_silicon)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxSiliconTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

   _seed_track_map  = findNode::getClass<TpcSeedTrackMap>(topNode, _tpcseed_track_map_name);
  if(!_seed_track_map)
    {
      std::cout << "Creating node TpcSeedTrackMap" << std::endl;

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
      PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
      
      /// Check that it is there
      if (!svtxNode)
	{
	  svtxNode = new PHCompositeNode("SVTX");
	  dstNode->addNode(svtxNode);
	}
      
      _seed_track_map = new TpcSeedTrackMapv1();
      PHIODataNode<PHObject> *node
	= new PHIODataNode<PHObject>(_seed_track_map, _tpcseed_track_map_name);
      svtxNode->addNode(node);
    }

 _corrected_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"CORRECTED_TRKR_CLUSTER");
  if(_corrected_cluster_map)
    {
      std::cout << " Found CORRECTED_TRKR_CLUSTER node " << std::endl;
    }
  
 _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
 if (!_cluster_map)
   {
     std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
     return Fun4AllReturnCodes::ABORTEVENT;
   }
 
 return Fun4AllReturnCodes::EVENT_OK;
}

void PHSiliconTpcTrackMatching::copySiliconClustersToCorrectedMap( )
{
  // loop over final track map, copy silicon clusters to corrected cluster map
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      SvtxTrack *track = phtrk_iter->second;

      // loop over associated clusters to get keys for silicon cluster
      for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
	   iter != track->end_cluster_keys();
	   ++iter)
	{
	  TrkrDefs::cluskey cluster_key = *iter;
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
	  if(trkrid == TrkrDefs::mvtxId || trkrid == TrkrDefs::inttId)
	    {
	      TrkrCluster *cluster =  _cluster_map->findCluster(cluster_key);	
	      TrkrCluster *newclus = _corrected_cluster_map->findOrAddCluster(cluster_key)->second;
	  
	      newclus->setSubSurfKey(cluster->getSubSurfKey());
	      newclus->setAdc(cluster->getAdc());
	      
	      newclus->setActsLocalError(0,0,cluster->getActsLocalError(0,0));
	      newclus->setActsLocalError(1,0,cluster->getActsLocalError(1,0));
	      newclus->setActsLocalError(0,1,cluster->getActsLocalError(0,1));
	      newclus->setActsLocalError(1,1,cluster->getActsLocalError(1,1));
	      
	      newclus->setLocalX(cluster->getLocalX());
	      newclus->setLocalY(cluster->getLocalY());
	    }
	}      
    }
}
  

