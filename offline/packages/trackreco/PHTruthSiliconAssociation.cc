#include "PHTruthSiliconAssociation.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>

#include "AssocInfoContainer.h"

using namespace std;

//____________________________________________________________________________..
PHTruthSiliconAssociation::PHTruthSiliconAssociation(const std::string &name):
 SubsysReco(name)
{
  //cout << "PHTruthSiliconAssociation::PHTruthSiliconAssociation(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHTruthSiliconAssociation::~PHTruthSiliconAssociation()
{

}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::Init(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);

  return ret;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1) 
    cout << "PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode) Processing Event" << endl;

  // Loop over all SvtxTracks from the CA seeder
  // These should contain all TPC clusters already

  const unsigned int original_track_map_lastkey = _track_map->end()->first;
  std::cout << " track map last key = " <<  original_track_map_lastkey << std::endl;

  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first >= original_track_map_lastkey)  break;

      _tracklet = phtrk_iter->second;
      if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet->get_phi()
	    << endl;
	}

      unsigned int vertexId = _tracklet->get_vertex_id();
      // if the vertex id from the seeder is nonsense, use vertex 0
      if(vertexId == UINT_MAX)
	vertexId = 0;
      _tracklet->set_vertex_id(vertexId);

      // set the track position to the vertex position
      const SvtxVertex *svtxVertex = _vertex_map->get(vertexId);
      
      _tracklet->set_x(svtxVertex->get_x());
      _tracklet->set_y(svtxVertex->get_y());
      _tracklet->set_z(svtxVertex->get_z());
      
      // identify the best truth track match(es) for this seed track       
      std::vector<PHG4Particle*> g4particle_vec = getG4PrimaryParticle(_tracklet);
      //std::cout << " g4particle_vec.size() " << g4particle_vec.size() << std::endl;

      if(g4particle_vec.size() < 1) continue;

      bool test_phi_matching = true;   // normally false
      if(test_phi_matching)
	{
	  // for getting the pT dependence of dphi  to eliminate the bias in phi from PHTpcTracker
	  if(g4particle_vec.size() == 1)
	    {
	      // print out the eta and phi values for analysis in case where match is unique
	      
	      double si_phi = atan2(g4particle_vec[0]->get_py(), g4particle_vec[0]->get_px());
	      double si_eta = asinh(g4particle_vec[0]->get_pz() / sqrt( pow(g4particle_vec[0]->get_px(),2)+pow( g4particle_vec[0]->get_py(),2)));
	      double si_pt = sqrt(pow(g4particle_vec[0]->get_px(),2)+pow( g4particle_vec[0]->get_py(),2));
	      
	      double tpc_phi = atan2(_tracklet->get_py(), _tracklet->get_px());
	      double tpc_eta = _tracklet->get_eta();
	      
	      cout << "         Try: " << "  pt " << si_pt << " tpc_phi " << tpc_phi << " si_phi " <<  si_phi 
		   << " tpc_eta " << tpc_eta << " si_eta " << si_eta << endl;
	      
	    }
	}
      
      // make copies of the original track for later use
      std::vector<SvtxTrack*> extraTrack; 
      for(unsigned int ig4=0;ig4 < g4particle_vec.size()-1; ++ig4)
	{      
	  SvtxTrack *newTrack = new SvtxTrack_v1();
	  // Not the first g4particle in the list, we need to add a new copy of the track to the track map and add the silicon clusters to that
	  const unsigned int lastTrackKey = _track_map->end()->first + ig4;
	  //std::cout << "   extra track key " << lastTrackKey << std::endl;
	 	  
	  newTrack->set_id(lastTrackKey);
	  newTrack->set_vertex_id(vertexId);

	  newTrack->set_charge(_tracklet->get_charge());
	  newTrack->set_px(_tracklet->get_px());
	  newTrack->set_py(_tracklet->get_py());
	  newTrack->set_pz(_tracklet->get_pz());
	  newTrack->set_x(_tracklet->get_x());
	  newTrack->set_y(_tracklet->get_y());
	  newTrack->set_z(_tracklet->get_z());
	  for(int i = 0; i < 6; ++i)
	    {
	      for(int j = 0; j < 6; ++j)
		{
		  newTrack->set_error(i,j, _tracklet->get_error(i,j));
		}
	    }

	  // loop over associated clusters to get hits for original TPC track, copy to new track
	  for (SvtxTrack::ConstClusterKeyIter iter = _tracklet->begin_cluster_keys();
	       iter != _tracklet->end_cluster_keys();
	       ++iter)
	    {
	      TrkrDefs::cluskey cluster_key = *iter;
	      newTrack->insert_cluster_key(cluster_key);
	    }
	  
	  extraTrack.push_back(newTrack);
	  //std::cout << "  added new copy of track " << lastTrackKey << " g4particle_vec.zize " << g4particle_vec.size() << std::endl;
	  //extraTrack[ig4]->identify();
	}

      // we just add the silicon clusters to the original track for the first returned g4particle
      
      if (Verbosity() >= 1) 
	{
	  std::cout << "  original track:" << endl;
	  _tracklet->identify(); 
	}
      
      // identify the clusters that are associated with this g4particle
      PHG4Particle* g4particle = g4particle_vec[0];
      std::set<TrkrDefs::cluskey> clusters = getSiliconClustersFromParticle(g4particle);

      for (std::set<TrkrDefs::cluskey>::iterator jter = clusters.begin();
	   jter != clusters.end();
	   ++jter)
	{
	  TrkrDefs::cluskey cluster_key = *jter;
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);
	  unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
	  if (Verbosity() >= 1)  std::cout << "     found silicon cluster with key: " << cluster_key << " in layer " << layer << std::endl;
	  
	  // Identify the MVTX and INTT clusters and add them to the SvtxTrack cluster key list
	  if(trkrid == TrkrDefs::mvtxId || trkrid == TrkrDefs::inttId)
	    {
	      if (Verbosity() >= 1)  std::cout << "            cluster belongs to MVTX or INTT, add to track " << std::endl;
	      _tracklet->insert_cluster_key(cluster_key);
	      _assoc_container->SetClusterTrackAssoc(cluster_key, _tracklet->get_id());
	    }
	}
      
      if (Verbosity() >= 1)
	{
	  std::cout << " updated original track:" << std::endl;
	  _tracklet->identify(); 
	  std::cout << " new cluster keys size " << _tracklet->size_cluster_keys() << endl;
	  std::cout << "Done with original track " << phtrk_iter->first << std::endl;
	}
      
      // now add any extra copies of the track      
      if(g4particle_vec.size() > 1)
	{
	  for(unsigned int ig4=0;ig4 < g4particle_vec.size()-1; ++ ig4)
	    {      
	      PHG4Particle* g4particle = g4particle_vec[ig4];
	      
	      if (Verbosity() >= 1)
		std::cout << "   ig4  " << ig4 << " g4particleID " << g4particle->get_track_id() 
			  << " px " <<  g4particle->get_px() 
			  << " py " <<  g4particle->get_py() 
			  << " phi " << atan2(g4particle->get_py(), g4particle->get_px() )
			  << std::endl;
	      
	      // identify the clusters that are associated with this g4particle
	      std::set<TrkrDefs::cluskey> clusters = getSiliconClustersFromParticle(g4particle);
	      	
	      // Not the first g4particle in the list, we need to add a new copy of the track to the track map and add the silicon clusters to that
	      const unsigned int lastTrackKey = _track_map->end()->first;
	      //std::cout << "   lastTrackKey " << lastTrackKey << std::endl;
	      
	      extraTrack[ig4]->set_id(lastTrackKey);
	      
	      if (Verbosity() >= 1) 
		{
		  std::cout << "  original (copy) track:" << endl;
		  extraTrack[ig4]->identify(); 
		}
	      
	      _track_map->insert(extraTrack[ig4]);
	      
	      for (std::set<TrkrDefs::cluskey>::iterator jter = clusters.begin();
		   jter != clusters.end();
		   ++jter)
		{
		  TrkrDefs::cluskey cluster_key = *jter;
		  unsigned int layer = TrkrDefs::getLayer(cluster_key);
		  unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
		  if (Verbosity() >= 1)  std::cout << "     found cluster with key: " << cluster_key << " in layer " << layer << std::endl;
		  
		  // Identify the MVTX and INTT clusters and add them to the SvtxTrack cluster key list
		  if(trkrid == TrkrDefs::mvtxId || trkrid == TrkrDefs::inttId)
		    {
		      if (Verbosity() >= 1)  std::cout << "            cluster belongs to MVTX or INTT, add to track " << std::endl;
		      extraTrack[ig4]->insert_cluster_key(cluster_key);
		      _assoc_container->SetClusterTrackAssoc(cluster_key, extraTrack[ig4]->get_id());
		    }
		}
	      
	      if (Verbosity() >= 1)
		{
		  std::cout << " updated additional track:" << std::endl;
		  extraTrack[ig4]->identify(); 
		  std::cout << " new cluster keys size " << extraTrack[ig4]->size_cluster_keys() << endl;
		  std::cout << "Done with extra track " << extraTrack[ig4]->get_id()  << std::endl;
		}
	    }
	}
    }
  
  if (Verbosity() >= 1)
    cout << "PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::ResetEvent(PHCompositeNode *topNode)
{
  //cout << "PHTruthSiliconAssociation::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::EndRun(const int runnumber)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::Reset(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PHTruthSiliconAssociation::Print(const std::string &what) const
{
  //cout << "PHTruthSiliconAssociation::Print(const std::string &what) const Printing info for " << what << endl;
}

int  PHTruthSiliconAssociation::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  _track_map = findNode::getClass<SvtxTrackMap>(topNode,  "SvtxTrackMap");
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
  if (!_assoc_container)
  {
    cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!_hit_truth_map)
  {
    cerr << PHWHERE << " ERROR: Can't find TrkrHitTruthAssoc." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!_cluster_hit_map)
  {
    cerr << PHWHERE << " ERROR: Can't find TrkrClusterHitAssoc." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

   _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<PHG4Particle*> PHTruthSiliconAssociation::getG4PrimaryParticle(SvtxTrack *track)
{
  // Find the best g4particle match to the clusters associated with this reco track

  std::vector<PHG4Particle*> g4part_vec;
  std::set<int> pid;
  std::multimap<int, int> pid_count;
  int minfound = 200;  // require at least this many hits to be put on the list of contributing particles

  // loop over associated clusters to get hits
  for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
       iter != track->end_cluster_keys();
       ++iter)
    {
      TrkrDefs::cluskey cluster_key = *iter;

      // get all reco hits for this cluster
      TrkrClusterHitAssoc::ConstRange hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey
      for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
	{
	  TrkrDefs::hitkey hitkey = clushititer->second;
	  // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
	  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);	  
	  
	  // get all of the g4hits for this hitkey
	  std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
	  _hit_truth_map->getG4Hits(hitsetkey, hitkey, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
	  for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); htiter != temp_map.end(); ++htiter) 
	    {
	      // extract the g4 hit key 
	      PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	      PHG4Hit * g4hit = nullptr;
	      unsigned int trkrid = TrkrDefs::getTrkrId(hitsetkey);
	      switch( trkrid )
		{
		case TrkrDefs::tpcId: g4hit = _g4hits_tpc->findHit(g4hitkey); break;
		case TrkrDefs::inttId: g4hit = _g4hits_intt->findHit(g4hitkey); break;
		case TrkrDefs::mvtxId: g4hit = _g4hits_mvtx->findHit(g4hitkey); break;
		default: break;
		}
	      if( g4hit )
		{
		  // get the g4particle ID for this g4hit
		  pid.insert(g4hit->get_trkid());
		  // this is for counting the number of times this g4particle was associated
		  int cnt = 1;
		  pid_count.insert(std::make_pair(g4hit->get_trkid(), cnt));
		}
	    } // end loop over g4hits associated with hitsetkey and hitkey
	} // end loop over hits associated with cluskey  
    }  // end loop over clusters associated with this track


  // loop over the particle id's, and count the number for each one
  for( auto it = pid.begin(); it != pid.end(); ++it)
    {
      if(*it < 0) continue;   // ignore secondary particles

      int nfound = 0;
      std::pair<std::multimap<int, int>::iterator, std::multimap<int,int>::iterator> this_pid = pid_count.equal_range(*it);
      for(auto cnt_it = this_pid.first; cnt_it != this_pid.second; ++cnt_it)
	{
	  nfound++;	  
	}

      if(Verbosity() >= 1) std::cout << "    pid: " << *it << "  nfound " << nfound << std::endl;
      if(nfound > minfound)
	{  
	  g4part_vec.push_back(_truthinfo->GetParticle(*it));
	}
    }

  return g4part_vec;
  
}

std::set<TrkrDefs::cluskey> PHTruthSiliconAssociation::getSiliconClustersFromParticle(PHG4Particle* g4particle)
{
  // Find the reco clusters in the silicon layers from this g4particle

  std::set<TrkrDefs::cluskey> clusters;

  // loop over all the clusters
  TrkrClusterContainer::ConstRange all_clusters = _cluster_map->getClusters();
  for (TrkrClusterContainer::ConstIterator iter = all_clusters.first;
       iter != all_clusters.second;
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = iter->first;
    unsigned int layer = TrkrDefs::getLayer(cluster_key);

    if(layer > 6) continue;  // we need the silicon layers only

    // get all truth hits for this cluster
    TrkrClusterHitAssoc::ConstRange hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey
    for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
      {
	TrkrDefs::hitkey hitkey = clushititer->second;
	  // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
	  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);	  
	  
	  // get all of the g4hits for this hitkey
	  std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
	  _hit_truth_map->getG4Hits(hitsetkey, hitkey, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
	  for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); htiter != temp_map.end(); ++htiter) 
	    {
	      // extract the g4 hit key 
	      PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	      PHG4Hit * g4hit = nullptr;
	      unsigned int trkrid = TrkrDefs::getTrkrId(hitsetkey);
	      switch( trkrid )
		{
		case TrkrDefs::mvtxId: g4hit = _g4hits_mvtx->findHit(g4hitkey); break;
		case TrkrDefs::inttId: g4hit = _g4hits_intt->findHit(g4hitkey); break;
		default: break;
		}
	      if( g4hit )
		{
		  // get the g4particle for this g4hit
		  int trkid = g4hit->get_trkid();	
		  if(trkid == g4particle->get_track_id())
		    {
		      clusters.insert(cluster_key);
		    }		    
		}
	    } // end loop over g4hits associated with hitsetkey and hitkey
      } // end loop over hits associated with cluskey  
  } // end loop over all cluster keys in silicon layers
  
  return clusters;
  
}

