#include "PHTruthSiliconAssociation.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterContainerv3.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterHitAssocv2.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack_v3.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase/TpcSeedTrackMapv1.h>    
#include <trackbase/TrkrClusterCrossingAssoc.h> 

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

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
int PHTruthSiliconAssociation::Init(PHCompositeNode */*topNode*/)
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
int PHTruthSiliconAssociation::process_event(PHCompositeNode */*topNode*/)
{
  if (Verbosity() >= 1) 
    cout << "PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode) Processing Event" << endl;

  // Loop over all SvtxTracks from the CA seeder
  // These should contain all TPC clusters already

  const unsigned int original_track_map_lastkey = _track_map->empty() ? 0:std::prev(_track_map->end())->first;
  if(Verbosity() > 0) std::cout << " track map last key = " <<  original_track_map_lastkey << std::endl;

  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      // we may add tracks to the map, so we stop at the last original track
      if(phtrk_iter->first > original_track_map_lastkey)  break;

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

      // this one is always there
      _seed_track_map->addAssoc(_tracklet->get_id(), _tracklet->get_id() ) ;

      // identify the best truth track match(es) for this seed track       
      std::vector<PHG4Particle*> g4particle_vec = getG4PrimaryParticle(_tracklet);
      if(Verbosity() > 0)  std::cout << " g4particle_vec.size() " << g4particle_vec.size() << std::endl;

      if(g4particle_vec.size() < 1) continue;

      bool test_phi_matching = false;   // normally false
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
	  //	  float time = g4particle_vec[ig4]->

	  SvtxTrack *newTrack = new SvtxTrack_v3();
	  // Not the first g4particle in the list, we need to add a new copy of the track to the track map and add the silicon clusters to that
	  const unsigned int lastTrackKey = ( _track_map->empty() ? 0:std::prev(_track_map->end())->first ) + ig4;
	  //std::cout << "   extra track key " << lastTrackKey + 1 << std::endl;
	 	  
	  newTrack->set_id(lastTrackKey + 1);
	  _seed_track_map->addAssoc(_tracklet->get_id(), newTrack->get_id()) ;

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

	      const int vertexID = g4particle->get_vtx_id();
	      
	      if (Verbosity() >= 1)
		std::cout << "   ig4  " << ig4 << " g4particleID " << g4particle->get_track_id() 
			  << " px " <<  g4particle->get_px() 
			  << " py " <<  g4particle->get_py() 
			  << " phi " << atan2(g4particle->get_py(), g4particle->get_px() )
			  << std::endl;
	      
	      // identify the clusters that are associated with this g4particle
	      std::set<TrkrDefs::cluskey> clusters = getSiliconClustersFromParticle(g4particle);
	      	
	      // Not the first g4particle in the list, we need to add a new copy of the track to the track map and add the silicon clusters to that
	      const unsigned int lastTrackKey =  _track_map->empty() ? 0:std::prev(_track_map->end())->first;
	      //std::cout << "   lastTrackKey " << lastTrackKey + 1 << std::endl;
	      
	      extraTrack[ig4]->set_id(lastTrackKey + 1);
	      
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
		    }
		}

	      // set the combined track vertex to the truth vertex for the silicon truth track
	      double random = 0;//((double) rand() / (RAND_MAX)) * 0.05;
	      // make it negative sometimes
	      if(rand() % 2)
		random *= -1;
	      // assign the track position using the truth vertex for this track
	      auto g4vertex = _g4truth_container->GetVtx(vertexID);
	      extraTrack[ig4]->set_x(g4vertex->get_x() * (1 + random));
	      extraTrack[ig4]->set_y(g4vertex->get_y() * (1 + random));
	      extraTrack[ig4]->set_z(g4vertex->get_z() * (1 + random));

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

  // loop over all tracks and copy the silicon clusters to the corrected cluster map
  if(_corrected_cluster_map)
    copySiliconClustersToCorrectedMap();

  // loop over all tracks and get the bunch crossing from the INTT cluster - crossing map and add it to the track
  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      SvtxTrack *track = phtrk_iter->second;
      std::vector<short int > intt_crossings = getInttCrossings(track);
      if(intt_crossings.size() == 0) 
	{
	  if(Verbosity() > 1) std::cout << " Silicon track " << track->get_id() << " has no INTT clusters" << std::endl;
	  continue ;
	}

      short int crossing_keep = intt_crossings[0];
      bool keep_it = true;
      for(unsigned int ic=1; ic<intt_crossings.size(); ++ic)
	{	  
	  if(intt_crossings[ic] != crossing_keep)
	    {
	      if(Verbosity() > 1) 
		std::cout << " INTT crossings not all the same for track " << track->get_id() << " crossing_keep " 
			  << crossing_keep << " new crossing " << intt_crossings[ic] << "- dropping this match " << std::endl;
	      keep_it = false;	      
	    }
	}
      if(keep_it)
	{            
	  track->set_crossing(crossing_keep);
	  if(Verbosity() > 1) std::cout << "                    Combined track " << track->get_id()  << " bunch crossing " << crossing_keep  << std::endl;           
	}
    }
  

  if (Verbosity() >= 1)
    cout << "PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::ResetEvent(PHCompositeNode */*topNode*/)
{
  //cout << "PHTruthSiliconAssociation::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::Reset(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PHTruthSiliconAssociation::Print(const std::string &/*what*/) const
{
  //cout << "PHTruthSiliconAssociation::Print(const std::string &what) const Printing info for " << what << endl;
}

int  PHTruthSiliconAssociation::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_g4truth_container)
    {
      cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _cluster_crossing_map = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
  if (!_cluster_crossing_map)
  {
    cerr << PHWHERE << " ERROR: Can't find TRKR_CLUSTERCROSSINGASSOC " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _corrected_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"CORRECTED_TRKR_CLUSTER");
  if(_corrected_cluster_map)
    {
      std::cout << " Found CORRECTED_TRKR_CLUSTER node " << std::endl;
    }
  
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!_hitsets)
    {
      std::cout << PHWHERE << "No hitset container on node tree. Bailing."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode,  "SvtxTrackMap");
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << endl;
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
  
  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  
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
      //TrkrClusterHitAssoc::ConstRange 
      std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>
	hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey
      //for(TrkrClusterHitAssoc::ConstIterator clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
      for(std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
	    clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
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
	  g4part_vec.push_back(_g4truth_container->GetParticle(*it));
	}
    }

  return g4part_vec;
  
}

std::set<TrkrDefs::cluskey> PHTruthSiliconAssociation::getSiliconClustersFromParticle(PHG4Particle* g4particle)
{
  // Find the reco clusters in the silicon layers from this g4particle

  std::set<TrkrDefs::cluskey> clusters;

  // loop over all the clusters
  auto hitsetrange = _hitsets->getHitSets();
  for (auto hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr){
    auto range = _cluster_map->getClusters(hitsetitr->first);
    for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){
      TrkrDefs::cluskey cluster_key = clusIter->first;
      unsigned int layer = TrkrDefs::getLayer(cluster_key);
      
      if(layer > 6) continue;  // we need the silicon layers only
      
      // get all truth hits for this cluster
      //TrkrClusterHitAssoc::ConstRange 
      std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>
	hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey
      //for(TrkrClusterHitAssoc::ConstIterator 
      for(std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
	    clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
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
  }//end loop over hitsets
  return clusters;
  
}

void PHTruthSiliconAssociation::copySiliconClustersToCorrectedMap( )
{
  // loop over final track map, copy silicon clusters to corrected cluster map
  for( auto track_iter = _track_map->begin(); track_iter != _track_map->end(); ++track_iter )
  {
    SvtxTrack* track = track_iter->second;
    // loop over associated clusters to get keys for micromegas cluster
    for(auto iter = track->begin_cluster_keys(); iter != track->end_cluster_keys(); ++iter)
    {
      TrkrDefs::cluskey cluster_key = *iter;
      const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
      if(trkrid == TrkrDefs::mvtxId || trkrid == TrkrDefs::inttId)
      {
        // check if clusters has not been inserted already
        if( _corrected_cluster_map->findCluster( cluster_key ) ) continue;

        auto cluster = _cluster_map->findCluster(cluster_key);	
        if( !cluster ) continue;
        
        // create a new cluster and copy from source
        auto newclus = new TrkrClusterv3;
        newclus->CopyFrom( cluster );

        // insert in corrected map
        _corrected_cluster_map->addCluster(newclus);
      }
    }      
  }
}

std::vector<short int> PHTruthSiliconAssociation::getInttCrossings(SvtxTrack *si_track)
{
  std::vector<short int> intt_crossings;

  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  // loop over associated clusters to get keys for silicon cluster
  for (SvtxTrack::ConstClusterKeyIter iter = si_track->begin_cluster_keys();
       iter != si_track->end_cluster_keys();
       ++iter)
    {
      
      TrkrDefs::cluskey cluster_key = *iter;
      const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
      if(trkrid == TrkrDefs::inttId)
	{
	  // std::cout << "      INTT cluster key " << cluster_key << std::endl; 
	  
	  TrkrCluster *cluster =  _cluster_map->findCluster(cluster_key);	
	  if( !cluster ) continue;	  
	  
	  unsigned int layer = TrkrDefs::getLayer(cluster_key);
	  
	  // get the bunch crossings for all hits in this cluster
	  auto crossings = _cluster_crossing_map->getCrossings(cluster_key);
	  for(auto iter = crossings.first; iter != crossings.second; ++iter)
	    {
	      std::cout << "      si Track " << si_track->get_id() << " cluster " << iter->first << " layer " << layer << " crossing " << iter->second  << std::endl;
	      intt_crossings.push_back(iter->second);
	    }
	}
    }
  //      std::cout << " intt_crossings size is " << intt_crossings.size() << std::endl;
  
  return intt_crossings;
}
