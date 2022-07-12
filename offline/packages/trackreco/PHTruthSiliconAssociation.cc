#include "PHTruthSiliconAssociation.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHRandomSeed.h>

/// Tracking includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/SvtxTrackSeed_v1.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase/TrkrClusterCrossingAssoc.h> 
#include <trackbase_historic/TrackSeed_FastSim_v1.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                            // for gsl_rng_alloc

using namespace std;

namespace
{ template< class T> inline constexpr T square( const T& x ) { return x*x; } }

//____________________________________________________________________________..
PHTruthSiliconAssociation::PHTruthSiliconAssociation(const std::string &name):
 SubsysReco(name)
{
  // initialize random generator
  const uint seed = PHRandomSeed();
  m_rng.reset( gsl_rng_alloc(gsl_rng_mt19937) );
  gsl_rng_set( m_rng.get(), seed );
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

  // Reset the silicon seed node
  _silicon_track_map->Reset();

  // Loop over all SeedTracks from the CA seeder
  // These should contain all TPC clusters already

  for (unsigned int phtrk_iter = 0;
       phtrk_iter < _tpc_track_map->size();
       ++phtrk_iter)
    {
      _tracklet = _tpc_track_map->get(phtrk_iter);

      if(!_tracklet)
	continue;
	
      if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter
	    << ": nhits: " << _tracklet-> size_cluster_keys()
	    << ": Total tracks: " << _tpc_track_map->size()
	    << endl;
	}

      // identify the best truth track match(es) for this seed track       
      std::vector<PHG4Particle*> g4particle_vec = getG4PrimaryParticle(_tracklet);
      if(Verbosity() > 0)  std::cout << " g4particle_vec.size() " << g4particle_vec.size() << std::endl;

      if(g4particle_vec.size() < 1) continue;

      if (Verbosity() >= 1) 
	{
	  std::cout << "  tpc seed track:" << endl;
	  _tracklet->identify(); 
	}

      for(unsigned int ig4=0;ig4 < g4particle_vec.size(); ++ig4)
	{      
	  // identify the clusters that are associated with this g4particle
	  PHG4Particle* g4particle = g4particle_vec[ig4];
	  std::set<TrkrDefs::cluskey> clusters = getSiliconClustersFromParticle(g4particle);
	  if(clusters.size() < 3) continue;

	  // Make a silicon track seed 
	  unsigned int silicon_seed_index = buildTrackSeed(clusters, g4particle, _silicon_track_map);


	  // Add this entry to the SvtxTrackSeedContainer	  
	  auto seed = std::make_unique<SvtxTrackSeed_v1>();
	  seed->set_tpc_seed_index(phtrk_iter);
	  seed->set_silicon_seed_index(silicon_seed_index);
	  _svtx_seed_map->insert(seed.get());

	  unsigned int combined_seed_index = _svtx_seed_map->size()-1;

	  if (Verbosity() >= 1)
	    {
	      std::cout << " Created SvtxTrackSeed  " << combined_seed_index 
			<< " with tpcid " << _svtx_seed_map->get(combined_seed_index)->get_tpc_seed_index()  
			<< " and silicon ID " << _svtx_seed_map->get(combined_seed_index)->get_silicon_seed_index() 
			<< std::endl;
	    }
	}
    }

  if(Verbosity() > 0)
    {
      // loop over the assembled tracks
      for (unsigned int phtrk_iter = 0;
	   phtrk_iter < _svtx_seed_map->size();
	   ++phtrk_iter)
	{
	  auto seed = _svtx_seed_map->get(phtrk_iter);
	  if(!seed) continue;
	  
	  auto tpc_index =  seed->get_tpc_seed_index();
	  auto silicon_index =  seed->get_silicon_seed_index();
	  
	  std::cout << "SvtxSeedTrack: " << phtrk_iter 
		    << " tpc_index " <<  tpc_index 
		    << " silicon_index " << silicon_index 
		    << std::endl;
	  
	  std::cout << " ----------  silicon tracklet " << silicon_index << std::endl;
	  auto silicon_tracklet = _silicon_track_map->get(silicon_index);
	  if(!silicon_tracklet) continue;
	  silicon_tracklet->identify();

	  std::cout << " ---------- tpc tracklet " << tpc_index << std::endl;
	  auto tpc_tracklet = _tpc_track_map->get(tpc_index);
	  if(!tpc_tracklet) continue;
	  tpc_tracklet->identify();
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

/*
void PHTruthSiliconAssociation::makeSvtxSeedMap()
{
  /// The track fitter expects a full svtxtrackseed container, so just 
  /// convert the tpc+silicon track seeds created into svtxtrack seeds
  
  for(unsigned int iter = 0; iter < _track_map->size(); ++iter)
    {
      auto seed = std::make_unique<SvtxTrackSeed_v1>();
      seed->set_tpc_seed_index(iter);
      _svtx_seed_map->insert(seed.get());
    }

}
*/

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
  /*
  _corrected_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode,"CORRECTED_TRKR_CLUSTER");
  if(_corrected_cluster_map)
    {
      std::cout << " Found CORRECTED_TRKR_CLUSTER node " << std::endl;
    }
  */
  
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
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
      
  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!_tgeometry)
    {
      std::cerr << PHWHERE << "ERROR: can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  /// Get the DST Node
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  
  /// Check that it is there
  if (!dstNode)
    {
      std::cerr << "DST Node missing, quitting" << std::endl;
      throw std::runtime_error("failed to find DST node in PHTruthSiliconAssociation::createNodes");
    }
  
  /// Get the tracking subnode under DST
  PHNodeIterator dstIter(dstNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));
  
  /// Check that it is there
  if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }

  _silicon_track_map = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
   if(!_silicon_track_map)
    {
      _silicon_track_map = new TrackSeedContainer_v1;
      PHIODataNode<PHObject>* si_tracks_node = 
	new PHIODataNode<PHObject>(_silicon_track_map, "SiliconTrackSeedContainer", "PHObject");
      svtxNode->addNode(si_tracks_node);
    }

  _tpc_track_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_tpc_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find TpcTrackSeedContainer: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  _svtx_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  if(!_svtx_seed_map)
    {
      _svtx_seed_map = new TrackSeedContainer_v1;
      PHIODataNode<PHObject> *node2 = new PHIODataNode<PHObject>(_svtx_seed_map, "SvtxTrackSeedContainer","PHObject");
      svtxNode->addNode(node2);
    }

  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<PHG4Particle*> PHTruthSiliconAssociation::getG4PrimaryParticle(TrackSeed *track)
{
  // Find the best g4particle match to the clusters associated with this reco track

  std::vector<PHG4Particle*> g4part_vec;
  std::set<int> pid;
  std::multimap<int, int> pid_count;
  int minfound = 200;  // require at least this many hits to be put on the list of contributing particles

  // loop over associated clusters to get hits
  for (auto iter = track->begin_cluster_keys();
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
  for(const auto& hitsetkey:_cluster_map->getHitSetKeys())
  {
    const auto layer = TrkrDefs::getLayer(hitsetkey);
    if(layer > 6) continue;  // we need the silicon layers only
    
    auto range = _cluster_map->getClusters(hitsetkey);
    for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){
      TrkrDefs::cluskey cluster_key = clusIter->first;
      
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

std::set<short int> PHTruthSiliconAssociation::getInttCrossings(TrackSeed *si_track) const
{
  std::set<short int> intt_crossings;

  // If the Si track contains an INTT hit, use it to get the bunch crossing offset
  // loop over associated clusters to get keys for silicon cluster
  for (auto iter = si_track->begin_cluster_keys();
       iter != si_track->end_cluster_keys();
       ++iter)
    {
      
      TrkrDefs::cluskey cluster_key = *iter;
      const unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
      if(trkrid == TrkrDefs::inttId)
	{
	  // std::cout << "      INTT cluster key " << cluster_key << std::endl; 
	  
	  const unsigned int layer = TrkrDefs::getLayer(cluster_key);
	  
	  // get the bunch crossings for all hits in this cluster
	  auto crossings = _cluster_crossing_map->getCrossings(cluster_key);
	  for(auto iter = crossings.first; iter != crossings.second; ++iter)
	    {
        const auto& [key, crossing] = *iter;
        if( Verbosity() )
        { std::cout << "PHTruthSiliconAssociation::getInttCrossings - si Track cluster " << key << " layer " << layer << " crossing " << crossing  << std::endl; }
	      intt_crossings.insert(crossing);
	    }
	}
    }
  //      std::cout << " intt_crossings size is " << intt_crossings.size() << std::endl;
  
  return intt_crossings;
}

unsigned int PHTruthSiliconAssociation::buildTrackSeed(std::set<TrkrDefs::cluskey> clusters, PHG4Particle *g4particle, TrackSeedContainer* container)
{
  auto track = std::make_unique<TrackSeed_FastSim_v1>();
  bool silicon = false;
  for (const auto& cluskey : clusters){
    if( TrkrDefs::getTrkrId(cluskey) == TrkrDefs::TrkrId::mvtxId || 
	TrkrDefs::getTrkrId(cluskey) == TrkrDefs::TrkrId::inttId)
      { silicon = true; }
    track->insert_cluster_key(cluskey);
  }
  
  const auto particle = TDatabasePDG::Instance()->GetParticle(g4particle->get_pid());
  int charge = 1;
  if(particle) 
    { 
      if(particle->Charge() < 0)
	{ charge = -1; }
    }

  auto random1 = gsl_ran_flat(m_rng.get(), 0.95, 1.05);  
  float px = g4particle->get_px() * random1;
  float py = g4particle->get_py() * random1;
  float pz = g4particle->get_pz() * random1;

  const auto g4vertex = _g4truth_container->GetVtx(g4particle->get_vtx_id());
  auto random2 = gsl_ran_flat(m_rng.get(), -0.002, 0.002);
  float x = g4vertex->get_x() + random2; 
  float y = g4vertex->get_y() + random2;
  float z = g4vertex->get_z() + random2;

  float pt = sqrt(px*px+py*py);
  float phi = atan2(py,px);
  float R = 100 * pt / (0.3*1.4);
  float theta = atan2(pt,pz);
  if(theta < 0)
    { theta += M_PI; }
  if(theta > M_PI)
    { theta -= M_PI; }
  
  float eta = -log(tan(theta/2.));

  // We have two equations, phi = atan2(-(X0-x),y-Y0) and 
  //R^2 = (x-X0)^2 + (y-Y0)^2. Solve for X0 and Y0 knowing R and phi
  float tanphisq = square(tan(phi));
  float a = tanphisq + 1;
  float b =-2*y*(tanphisq+1);
  float c = (tanphisq+1)*square(y)-square(R);
  
  float Y0_1 = (-b + sqrt(square(b)-4*a*c)) / (2.*a);
  float Y0_2 = (-b - sqrt(square(b)-4*a*c)) / (2.*a);
  float X0_1 = sqrt(pow(R, 2) - pow(Y0_1 - y, 2)) + x;
  float X0_2 = -sqrt(pow(R, 2) - pow(Y0_2 - y, 2)) + x;
  track->set_X0(X0_1);
  track->set_Y0(Y0_1);
  track->set_qOverR(charge / R);
  track->set_slope(1. / tan(theta));
  track->set_Z0(z);
  
  /// Need to find the right one for the bend angle
  
  float newphi = track->get_phi(_cluster_map, _tgeometry);
  /// We have to pick the right one based on the bend angle, so iterate
  /// through until you find the closest phi match
  if( fabs(newphi-phi) > 0.03)
    {
      track->set_X0(X0_2);
      newphi = track->get_phi(_cluster_map, _tgeometry);
  
      if( fabs(newphi-phi) > 0.03)
	{
	  track->set_Y0(Y0_2);
	  newphi = track->get_phi(_cluster_map, _tgeometry);

	  if( fabs(newphi-phi) > 0.03)
	    {
	      track->set_X0(X0_1);
	      newphi = track->get_phi(_cluster_map, _tgeometry);
	    }
	}
    }
  
  if(Verbosity() > 2)
    {
      std::cout << "Charge is " << charge << std::endl;
      std::cout << "truth/reco px " << px << ", " << track->get_px(_cluster_map, _tgeometry) << std::endl;
      std::cout << "truth/reco py " << py << ", " << track->get_py(_cluster_map, _tgeometry) << std::endl;
      std::cout << "truth/reco pz " << pz << ", " << track->get_pz() << std::endl;
      std::cout << "truth/reco pt " << pt << ", " << track->get_pt() << std::endl;
      std::cout << "truth/reco phi " << phi << ", " << track->get_phi(_cluster_map, _tgeometry) << std::endl;
      std::cout << "truth/reco eta " << eta << ", " << track->get_eta() << std::endl;
      std::cout << "truth/reco x " << x << ", " << track->get_x() << std::endl;
      std::cout << "truth/reco y " << y << ", " << track->get_y() << std::endl;
      std::cout << "truth/reco z " << z << ", " << track->get_z() << std::endl;
	
    }
  
  // set intt crossing
  if(silicon)
    {
      // silicon tracklet
      /* inspired from PHtruthSiliconAssociation */
      const auto intt_crossings = getInttCrossings(track.get());
      if(intt_crossings.empty()) 
	{
	  if(Verbosity() > 1)  std::cout << "PHTruthTrackSeeding::Process - Silicon track " << container->size() - 1 << " has no INTT clusters" << std::endl;
	  return (container->size() -1);
	} else if( intt_crossings.size() > 1 ) {
	if(Verbosity() > 1) 
	  { std::cout << "PHTruthTrackSeeding::Process - INTT crossings not all the same for track " << container->size() - 1 << " crossing_keep - dropping this match " << std::endl; }
	
      } else {
	const auto& crossing = *intt_crossings.begin();
	track->set_crossing(crossing);
	if(Verbosity() > 1)
	  std::cout << "PHTruthTrackSeeding::Process - Combined track " << container->size() - 1  << " bunch crossing " << crossing << std::endl;           
      }
    }  // end if _min_layer
  else
    {
      // no INTT layers, crossing is unknown
      track->set_crossing(SHRT_MAX);	
    }

  container->insert(track.get());

  return (container->size() - 1);
}
