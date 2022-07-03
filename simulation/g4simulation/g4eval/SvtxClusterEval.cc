#include "SvtxClusterEval.h"

#include "SvtxHitEval.h"
#include "SvtxTruthEval.h"

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase_historic/ActsTransformations.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/getClass.h>
#include <phool/PHTimer.h>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>                          // for operator<<, basic_ostream
#include <map>
#include <set>
#include <TVector3.h>

using namespace std;

SvtxClusterEval::SvtxClusterEval(PHCompositeNode* topNode)
  : _hiteval(topNode)
{
  get_node_pointers(topNode);
}

SvtxClusterEval::~SvtxClusterEval()
{
  if (_verbosity > 0)
    {
      if ((_errors > 0) || (_verbosity > 1))
	{
	  cout << "SvtxClusterEval::~SvtxClusterEval() - Error Count: " << _errors << endl;
	}
    }
}

void SvtxClusterEval::next_event(PHCompositeNode* topNode)
{
  _cache_all_truth_hits.clear();
  _cache_all_truth_clusters.clear();
  _cache_max_truth_hit_by_energy.clear();
  _cache_max_truth_cluster_by_energy.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_energy.clear();
  _cache_max_truth_particle_by_cluster_energy.clear();
  _cache_all_clusters_from_particle.clear();
  _cache_all_clusters_from_g4hit.clear();
  _cache_best_cluster_from_g4hit.clear();
  _cache_get_energy_contribution_g4particle.clear();
  _cache_get_energy_contribution_g4hit.clear();
  _cache_best_cluster_from_gtrackid_layer.clear();
  _clusters_per_layer.clear();
  //  _g4hits_per_layer.clear();
  _hiteval.next_event(topNode);
  
  get_node_pointers(topNode);
} 

std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > SvtxClusterEval::all_truth_clusters(TrkrDefs::cluskey cluster_key)
{
  if (_do_cache)
  {
    const auto iter = _cache_all_truth_clusters.find(cluster_key);
    if (iter != _cache_all_truth_clusters.end())
    {
      return iter->second;
    }
  }

  std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > truth_clusters;

  unsigned int cluster_layer = TrkrDefs::getLayer(cluster_key);
  std::set<PHG4Particle*> particles = all_truth_particles(cluster_key);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
    {
      PHG4Particle* particle = *iter;
      
      for( const auto& [ckey,cluster]:get_truth_eval()->all_truth_clusters(particle) )
      {
        if( TrkrDefs::getLayer(ckey) == cluster_layer )
        { truth_clusters.insert(std::make_pair(ckey, cluster)); }
      }
    }
  
    
  return _cache_all_truth_clusters.insert(std::make_pair(cluster_key, truth_clusters)).first->second;
}

std::pair<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> SvtxClusterEval::max_truth_cluster_by_energy(TrkrDefs::cluskey cluster_key)
{
  if (_do_cache)
  {
    const auto iter = _cache_max_truth_cluster_by_energy.find(cluster_key);
    if (iter != _cache_max_truth_cluster_by_energy.end())
    {
      return iter->second;
    }
  }

  unsigned int cluster_layer = TrkrDefs::getLayer(cluster_key);

  PHG4Particle* max_particle = max_truth_particle_by_cluster_energy(cluster_key);
  if(!max_particle) return std::make_pair( 0, nullptr );

  if(_verbosity > 0) 
    cout << "         max truth particle by cluster energy has  trackID  " << max_particle->get_track_id() << endl;      

  TrkrCluster* reco_cluster = _clustermap->findCluster(cluster_key);

  auto global = getGlobalPosition(cluster_key, reco_cluster);
  double reco_x = global(0);
  double reco_y = global(1);
  double reco_z = global(2);
  double r = sqrt(reco_x*reco_x + reco_y*reco_y);
  //double reco_rphi = r*fast_approx_atan2(reco_y, reco_x);
  double reco_rphi = r*atan2(reco_y, reco_x);
  
  const std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > gclusters = get_truth_eval()->all_truth_clusters(max_particle);
  for( const auto& [ckey,candidate_truth_cluster] : gclusters )
  {
    if(TrkrDefs::getLayer(ckey) != cluster_layer)  continue;

	  double gx = candidate_truth_cluster->getX();
	  double gy = candidate_truth_cluster->getY();
	  double gz = candidate_truth_cluster->getZ();
	  double gr = sqrt(gx*gx+gy*gy);
	  double grphi = gr*atan2(gy, gx);
	  //double grphi = gr*fast_approx_atan2(gy, gx);
	  
	  // Find the difference in position from the reco cluster
	  double dz = reco_z - gz;
	  double drphi = reco_rphi - grphi;
	  
	  // approximate 4 sigmas cut
	  if(cluster_layer > 6 && cluster_layer < 23)
	    {
	      if(fabs(drphi) < 4.0 * sig_tpc_rphi_inner &&
		 fabs(dz) < 4.0 * sig_tpc_z)
		{
		  return std::make_pair(ckey,candidate_truth_cluster);
		}
	    }
	  if(cluster_layer > 22 && cluster_layer < 39)
	    {
	      if(fabs(drphi) < 4.0 * sig_tpc_rphi_mid &&
		 fabs(dz) < 4.0 * sig_tpc_z)
		{
		  return std::make_pair(ckey,candidate_truth_cluster);
		}
	    }
	  if(cluster_layer > 38 && cluster_layer < 55)
	    {
	      if(fabs(drphi) < 4.0 * sig_tpc_rphi_outer &&
		 fabs(dz) < 4.0 * sig_tpc_z)
		{
		  return std::make_pair(ckey,candidate_truth_cluster);
		}
	    }
	  else if(cluster_layer < 3)
	    {
	      if(fabs(drphi) < 4.0 * sig_mvtx_rphi &&
		 fabs(dz) < 4.0 * sig_mvtx_z)
		{
		  return std::make_pair(ckey,candidate_truth_cluster);
		}
	    }
	  else if(cluster_layer == 55)
	    {
	      if(fabs(drphi) < 4.0 * sig_mms_rphi_55)
		{
		  return std::make_pair(ckey,candidate_truth_cluster);
		}
	    }
	  else if(cluster_layer == 56)
	    {
	      if(fabs(dz) < 4.0 * sig_mms_z_56)
		{
		  return std::make_pair(ckey,candidate_truth_cluster);
		}
	    }
	  else
	    {
	      if(fabs(drphi) < 4.0 * sig_intt_rphi &&
		 fabs(dz) < range_intt_z)
		{
		  return std::make_pair(ckey,candidate_truth_cluster);
		}
	    }
    }
  
		  return std::make_pair(0, nullptr);
}

std::pair<TrkrDefs::cluskey, TrkrCluster*> SvtxClusterEval::reco_cluster_from_truth_cluster(TrkrDefs::cluskey ckey, std::shared_ptr<TrkrCluster> gclus)
{
  if (_do_cache)
  {
    /* this does not work. Cache is not filled in the code below, so always remains empty */
    const auto iter = _cache_reco_cluster_from_truth_cluster.find(gclus);
    if (iter != _cache_reco_cluster_from_truth_cluster.end())
    {
      return iter->second;
    }
  }

  double gx = gclus->getX();
  double gy = gclus->getY();
  double gz = gclus->getZ();
  double gr = sqrt(gx*gx+gy*gy);
  double grphi = gr*atan2(gy, gx);
  //double grphi = gr*fast_approx_atan2(gy, gx);

  unsigned int truth_layer = TrkrDefs::getLayer(ckey);

  std::set<TrkrDefs::cluskey> reco_cluskeys;
  std::set<PHG4Hit*> contributing_hits =  get_truth_eval()->get_truth_hits_from_truth_cluster(ckey);
  for(auto it = contributing_hits.begin(); it != contributing_hits.end(); ++it)
    {
      PHG4Hit* cont_g4hit = *it;            
      std::set<TrkrDefs::cluskey> cluskeys = all_clusters_from(cont_g4hit);  // this returns clusters from this hit in any layer using TrkrAssoc maps
      
      if(_verbosity > 0)
	cout << "       contributing g4hitID " << cont_g4hit->get_hit_id() << " g4trackID " << cont_g4hit->get_trkid() << endl;
      
      for (std::set<TrkrDefs::cluskey>::iterator iter = cluskeys.begin();
	   iter != cluskeys.end();
	   ++iter)
	{
	  unsigned int clus_layer = TrkrDefs::getLayer(*iter);
	  if(clus_layer != truth_layer)  continue;
	  
	  reco_cluskeys.insert(*iter);
	}
    }

  unsigned int nreco = reco_cluskeys.size();
  if(nreco > 0)
    {
      // Find a matching reco cluster with position inside 4 sigmas, and replace reco_cluskey
      for( const auto& this_ckey:reco_cluskeys )
	{
	  // get the cluster
	  TrkrCluster* this_cluster = _clustermap->findCluster(this_ckey);
	  auto global = getGlobalPosition(this_ckey,this_cluster);
	  double this_x = global(0);
	  double this_y = global(1);
	  double this_z = global(2);
	  double this_rphi = gr*atan2(this_y, this_x);
	  //double this_rphi = gr*fast_approx_atan2(this_y, this_x);
	  
	  // Find the difference in position from the g4cluster
	  double dz = this_z - gz;
	  double drphi = this_rphi - grphi;

	  // approximate 4 sigmas cut
	  if(truth_layer > 6 && truth_layer < 23)
	    {
	      if(fabs(drphi) < 4.0 * sig_tpc_rphi_inner &&
		 fabs(dz) < 4.0 * sig_tpc_z)
		{
		  return std::make_pair( this_ckey, this_cluster);
		}
	    }
	  if(truth_layer > 22 && truth_layer < 39)
	    {
	      if(fabs(drphi) < 4.0 * sig_tpc_rphi_mid &&
		 fabs(dz) < 4.0 * sig_tpc_z)
		{
		  return std::make_pair( this_ckey, this_cluster);
		}
	    }
	  if(truth_layer > 38 && truth_layer < 55)
	    {
	      if(fabs(drphi) < 4.0 * sig_tpc_rphi_outer &&
		 fabs(dz) < 4.0 * sig_tpc_z)
		{
		  return std::make_pair( this_ckey, this_cluster);
		}
	    }
	  else if(truth_layer < 3)
	    {
	      if(fabs(drphi) < 4.0 * sig_mvtx_rphi &&
		 fabs(dz) < 4.0 * sig_mvtx_z)
		{
		  return std::make_pair( this_ckey, this_cluster);
		}
	    }
	  else if(truth_layer == 55)
	    {
	      if(fabs(drphi) < 4.0 * sig_mms_rphi_55)
		{
		  return std::make_pair( this_ckey, this_cluster);
		}
	    }
	  else if(truth_layer == 56)
	    {
	      if(fabs(dz) < 4.0 * sig_mms_z_56)
		{
		  return std::make_pair( this_ckey, this_cluster);
		}
	    }
	  else
	    {
	      if(fabs(drphi) < 4.0 * sig_intt_rphi &&
		 fabs(dz) < range_intt_z)
		{
		  return std::make_pair( this_ckey, this_cluster);
		}
	    }
	}
    } 
      
  return std::make_pair(0, nullptr);      
}

std::set<PHG4Hit*> SvtxClusterEval::all_truth_hits(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return std::set<PHG4Hit*>();
    }
  
  if (_do_cache)
    {
      std::map<TrkrDefs::cluskey, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(cluster_key);
      if (iter != _cache_all_truth_hits.end())
	{
	  return iter->second;
	}
    }
  
  std::set<PHG4Hit*> truth_hits;

  // get all truth hits for this cluster
  //_cluster_hit_map->identify();
  std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator> 
    hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey

  for(std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
	clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
    {
      TrkrDefs::hitkey hitkey = clushititer->second;
      // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey
      TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);	  
      
      // get all of the g4hits for this hitkey
      std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
      _hit_truth_map->getG4Hits(hitsetkey, hitkey, temp_map); 	  
      // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only

      for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); htiter != temp_map.end(); ++htiter) 
	{
	  
	  // extract the g4 hit key here and add the hits to the set
	  PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	  PHG4Hit * g4hit = nullptr;
	  unsigned int trkrid = TrkrDefs::getTrkrId(hitsetkey);
	  switch( trkrid )
	    {
	    case TrkrDefs::tpcId: g4hit = _g4hits_tpc->findHit(g4hitkey); break;
	    case TrkrDefs::inttId: g4hit = _g4hits_intt->findHit(g4hitkey); break;
	    case TrkrDefs::mvtxId: g4hit = _g4hits_mvtx->findHit(g4hitkey); break;
	    case TrkrDefs::micromegasId: g4hit = _g4hits_mms->findHit(g4hitkey); break;
	    default: break;
	    }
	  if( g4hit ) truth_hits.insert(g4hit);	      
	} // end loop over g4hits associated with hitsetkey and hitkey
    } // end loop over hits associated with cluskey  
  
  if (_do_cache) _cache_all_truth_hits.insert(make_pair(cluster_key, truth_hits));

  return truth_hits;
}

PHG4Hit* SvtxClusterEval::all_truth_hits_by_nhit(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return nullptr;
    }
  
//  if (_strict)
//    {
//      assert(cluster_key);
//    }
//  else if (!cluster_key)
//    {
//      ++_errors;
//      return std::set<PHG4Hit*>();
//    }
  /*
  if (_do_cache)
    {
      std::map<TrkrDefs::cluskey, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(cluster_key);
      if (iter != _cache_all_truth_hits.end())
	{
	  return iter->second;
	}
    }
  */
  TrkrCluster* cluster = _clustermap->findCluster(cluster_key);
  auto glob = getGlobalPosition(cluster_key, cluster);
  TVector3 cvec(glob(0), glob(1), glob(2));
  unsigned int layer = TrkrDefs::getLayer(cluster_key);
  std::set<PHG4Hit*> truth_hits;

  std::multimap<PHG4HitDefs::keytype,TrkrDefs::hitkey> g4keyperhit;
  std::vector<PHG4HitDefs::keytype> g4hitkeys;
  // get all truth hits for this cluster
  //_cluster_hit_map->identify();
  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);	  

  std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator> 
    hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey
  for(std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
	clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
    {
      TrkrDefs::hitkey hitkey = clushititer->second;
      // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey

      // get all of the g4hits for this hitkey
      std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
      _hit_truth_map->getG4Hits(hitsetkey, hitkey, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
      for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); htiter != temp_map.end(); ++htiter) 
	{
	  // extract the g4 hit key here and add the hits to the set
	  PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	  if(_verbosity > 2)
	    cout << " g4key:  " <<  g4hitkey << " layer: " << layer << endl;
	  TrkrDefs::hitkey hitkey =  htiter->second.first;
	  /*	  if(layer>=7){
	    PHG4Hit *match_g4hit = _g4hits_tpc->findHit(g4hitkey);
	    if(layer != match_g4hit->get_layer() ) continue;
	    }
	  */
	  g4keyperhit.insert(std::pair<PHG4HitDefs::keytype,TrkrDefs::hitkey>(g4hitkey,hitkey));
	  std::vector<PHG4HitDefs::keytype>::iterator itg4keys = find(g4hitkeys.begin(),g4hitkeys.end(),g4hitkey);
	  if(itg4keys==g4hitkeys.end()) g4hitkeys.push_back(g4hitkey);
	} // end loop over g4hits associated with hitsetkey and hitkey
    } // end loop over hits associated with cluskey  

  //  if (_do_cache) _cache_all_truth_hits.insert(make_pair(cluster_key, truth_hits));
  PHG4HitDefs::keytype max_key = 0;
  unsigned int n_max = 0;

  if(g4hitkeys.size()==1 ){
    std::vector<PHG4HitDefs::keytype>::iterator it = g4hitkeys.begin();
    max_key = *it;
  }else{
    for(std::vector<PHG4HitDefs::keytype>::iterator it = g4hitkeys.begin(); it != g4hitkeys.end(); ++it){
      unsigned int ng4hit = g4keyperhit.count(*it);
      PHG4Hit * this_g4hit = _g4hits_tpc->findHit(*it);
      
      if(layer >= 7 ){ //in tpc
	if(this_g4hit!=NULL){
	  unsigned int glayer = this_g4hit->get_layer();
	  if(layer != glayer) continue;
	  
	  TVector3 vec(this_g4hit->get_avg_x(), this_g4hit->get_avg_y(), this_g4hit->get_avg_z());
	  //cout << "layer: " << layer << " (" << glayer << ") " << " gtrackID: " << this_g4hit->get_trkid() << " novlp: " << ng4hit << " phi: " << vec.Phi() << " z: " << this_g4hit->get_avg_z() << " r: " << vec.Perp() << " keyg4: " << *it << endl; //<< " keyrec: "<< *it.second << endl;
	}
	/*else{
	  cout << "g4hit == NULL " << endl; 
	}
	*/
      }
      if(ng4hit>n_max){
	max_key = *it;
	n_max = ng4hit;
      }
      
    }
  }

  PHG4Hit * g4hit = nullptr;
  unsigned int trkrid = TrkrDefs::getTrkrId(hitsetkey);
  switch( trkrid )
    {
    case TrkrDefs::tpcId: g4hit        = _g4hits_tpc->findHit(max_key); break;
    case TrkrDefs::inttId: g4hit       = _g4hits_intt->findHit(max_key); break;
    case TrkrDefs::mvtxId: g4hit       = _g4hits_mvtx->findHit(max_key); break;
    case TrkrDefs::micromegasId: g4hit = _g4hits_mms->findHit(max_key); break;
    default: break;
    }
  if( g4hit ) truth_hits.insert(g4hit);	 

  return g4hit;
}

std::pair<int, int> SvtxClusterEval::gtrackid_and_layer_by_nhit(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return make_pair(0,0);
    }
  
//  if (_strict)
//    {
//      assert(cluster_key);
//    }
//  else if (!cluster_key)
//    {
//      ++_errors;
//      return std::set<PHG4Hit*>();
//    }
  /*
  if (_do_cache)
    {
      std::map<TrkrDefs::cluskey, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(cluster_key);
      if (iter != _cache_all_truth_hits.end())
	{
	  return iter->second;
	}
    }
  */

  std::pair<int, int> out_pair;
  out_pair.first = 0;
  out_pair.second = -1;

  TrkrCluster* cluster = _clustermap->findCluster(cluster_key);
  auto global = getGlobalPosition(cluster_key, cluster);
  TVector3 cvec(global(0), global(1), global(2));
  unsigned int layer = TrkrDefs::getLayer(cluster_key);

  std::multimap<PHG4HitDefs::keytype,TrkrDefs::hitkey> g4keyperhit;
  std::vector<PHG4HitDefs::keytype> g4hitkeys;
  // get all truth hits for this cluster
  //_cluster_hit_map->identify();
  TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);	  

  std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator> 
    hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey
  for(std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
	clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
    {
      TrkrDefs::hitkey hitkey = clushititer->second;
      // TrkrHitTruthAssoc uses a map with (hitsetkey, std::pair(hitkey, g4hitkey)) - get the hitsetkey from the cluskey

      // get all of the g4hits for this hitkey
      std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> > temp_map;    
      _hit_truth_map->getG4Hits(hitsetkey, hitkey, temp_map); 	  // returns pairs (hitsetkey, std::pair(hitkey, g4hitkey)) for this hitkey only
    
      for(std::multimap< TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> >::iterator htiter =  temp_map.begin(); htiter != temp_map.end(); ++htiter) 
	{
	  // extract the g4 hit key here and add the hits to the set
	  PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	  if(_verbosity > 2)
	    cout << " g4key:  " <<  g4hitkey << " layer: " << layer << endl;
	  TrkrDefs::hitkey hitkey =  htiter->second.first;
	  /*	  if(layer>=7){
	    PHG4Hit *match_g4hit = _g4hits_tpc->findHit(g4hitkey);
	    if(layer != match_g4hit->get_layer() ) continue;
	    }
	  */
	  g4keyperhit.insert(std::pair<PHG4HitDefs::keytype,TrkrDefs::hitkey>(g4hitkey,hitkey));
	  std::vector<PHG4HitDefs::keytype>::iterator itg4keys = find(g4hitkeys.begin(),g4hitkeys.end(),g4hitkey);
	  if(itg4keys==g4hitkeys.end()) g4hitkeys.push_back(g4hitkey);
	} // end loop over g4hits associated with hitsetkey and hitkey
    } // end loop over hits associated with cluskey  

  PHG4HitDefs::keytype max_key = 0;
  unsigned int n_max = 0;
  if(_verbosity > 2)
    cout << " n matches found: " << g4hitkeys.size() << " phi: " << cvec.Phi() << " z: " << cvec.Z() << " ckey: " << cluster_key << endl;

  if(g4hitkeys.size()==1 ){
    std::vector<PHG4HitDefs::keytype>::iterator it = g4hitkeys.begin();
    max_key = *it;
  }else{
    for(std::vector<PHG4HitDefs::keytype>::iterator it = g4hitkeys.begin(); it != g4hitkeys.end(); ++it){
      unsigned int ng4hit = g4keyperhit.count(*it);
      PHG4Hit * this_g4hit = _g4hits_tpc->findHit(*it);
      
      
      if(layer >= 7 ){ //in tpc
	if(this_g4hit!=NULL){
	  unsigned int glayer = this_g4hit->get_layer();
	  //  if(layer != glayer) continue;
	    
	  TVector3 vec(this_g4hit->get_avg_x(), this_g4hit->get_avg_y(), this_g4hit->get_avg_z());
	  if(_verbosity > 2)
	    cout << "layer: " << layer << " (" << glayer << ") " << " gtrackID: " << this_g4hit->get_trkid() << " novlp: " << ng4hit << " phi: " << vec.Phi() << " z: " << this_g4hit->get_avg_z() << " r: " << vec.Perp() << " keyg4: " << *it << " cz: " << cluster->getZ() << endl; //<< " keyrec: "<< *it.second << endl;
	}
      }
      if(ng4hit>n_max){
	max_key = *it;
	n_max = ng4hit;
      }
      
    }
  }
  if(_verbosity > 2)
    cout << "found in layer: " << layer << " n_max: " << n_max << " max_key: " << max_key << " ckey: " << cluster_key << endl;  
  if(max_key !=0){
    PHG4Hit * g4hit = nullptr;
    unsigned int trkrid = TrkrDefs::getTrkrId(hitsetkey);
    switch( trkrid )
    {
    case TrkrDefs::tpcId: g4hit = _g4hits_tpc->findHit(max_key); break;
    case TrkrDefs::inttId: g4hit = _g4hits_intt->findHit(max_key); break;
    case TrkrDefs::mvtxId: g4hit = _g4hits_mvtx->findHit(max_key); break;
    case TrkrDefs::micromegasId: g4hit = _g4hits_mms->findHit(max_key); break;
    default: break;
    }
    
    //check if we on a looper
    PHG4Particle* g4particle = _truthinfo->GetParticle(g4hit->get_trkid());

    PHG4VtxPoint* vtx = _truthinfo->GetVtx(g4particle->get_vtx_id());
    float vtx_z = vtx->get_z();
    float gpx = g4particle->get_px();
    float gpy = g4particle->get_py();
    float gpz = g4particle->get_pz();
    float gpeta = NAN;
    
    TVector3 gv(gpx, gpy, gpz);
    gpeta = gv.Eta();
    TVector3 this_vec( g4hit->get_avg_x() ,
		       g4hit->get_avg_y() ,
		       g4hit->get_avg_z() - vtx_z);
    double deta = TMath::Abs(gpeta - this_vec.Eta());
    
    int is_loop = 0;
    
    if(layer >= 7){
      //	    cout << " in tpc " << endl;
      if(deta>0.1) is_loop = 1;
    }
    
    out_pair.first = g4hit->get_trkid();
    if(!is_loop)
      out_pair.second = layer;
  }
  return out_pair;
}

PHG4Hit* SvtxClusterEval::max_truth_hit_by_energy(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return nullptr;
    }
  
  if (_do_cache)
    {
      std::map<TrkrDefs::cluskey, PHG4Hit*>::iterator iter =
        _cache_max_truth_hit_by_energy.find(cluster_key);
      if (iter != _cache_max_truth_hit_by_energy.end())
	{
	  return iter->second;
	}
    }
  
  std::set<PHG4Hit*> hits = all_truth_hits(cluster_key);
  PHG4Hit* max_hit = nullptr;
  float max_e = FLT_MAX * -1.0;
  for (std::set<PHG4Hit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter)
    {
      PHG4Hit* hit = *iter;
      if (hit->get_edep() > max_e)
	{
	  max_e = hit->get_edep();
	  max_hit = hit;
	}
    }
  
  if (_do_cache) _cache_max_truth_hit_by_energy.insert(make_pair(cluster_key, max_hit));
  
  return max_hit;
}

std::set<PHG4Particle*> SvtxClusterEval::all_truth_particles(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return std::set<PHG4Particle*>();
    }
  
  if (_do_cache)
    {
      std::map<TrkrDefs::cluskey, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(cluster_key);
      if (iter != _cache_all_truth_particles.end())
	{
	  return iter->second;
	}
    }
  
  std::set<PHG4Particle*> truth_particles;
  
  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster_key);
  
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter)
    {
      PHG4Hit* hit = *iter;
      PHG4Particle* particle = get_truth_eval()->get_particle(hit);
      //cout << "cluster key " << cluster_key << " has hit " << hit->get_hit_id() << " and has particle " << particle->get_track_id() << endl;

      if (_strict)
	{
	  assert(particle);
	}
      else if (!particle)
	{
	  ++_errors;
	  continue;
	}

      truth_particles.insert(particle);
    }
  
  if (_do_cache) _cache_all_truth_particles.insert(make_pair(cluster_key, truth_particles));
  
  return truth_particles;
}

PHG4Particle* SvtxClusterEval::max_truth_particle_by_cluster_energy(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return nullptr;
    }
   
  if (_do_cache)
    {
      std::map<TrkrDefs::cluskey, PHG4Particle*>::iterator iter =
        _cache_max_truth_particle_by_cluster_energy.find(cluster_key);
      if (iter != _cache_max_truth_particle_by_cluster_energy.end())
	{
	  return iter->second;
	}
    }

  unsigned int layer = TrkrDefs::getLayer(cluster_key);
  
  // loop over all particles associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(cluster_key);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
    {
      PHG4Particle* particle = *iter;
      std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > truth_clus = get_truth_eval()->all_truth_clusters(particle);
      for( const auto& [ckey, cluster]:truth_clus )
	{	
    if( TrkrDefs::getLayer(ckey) == layer)
	    {
	      float e = cluster->getError(0,0);
        if (e > max_e)
		{
		  max_e = e;
		  max_particle = particle;
		}
	    }
	}
    }
 
  if (_do_cache) _cache_max_truth_particle_by_cluster_energy.insert(make_pair(cluster_key, max_particle));
  
  return max_particle;
}

PHG4Particle* SvtxClusterEval::max_truth_particle_by_energy(TrkrDefs::cluskey cluster_key)
{
  // Note: this does not quite work correctly for the TPC - it assumes one g4hit per layer
  // use max_truth_particle_by_cluster_energy instead
 
  if (!has_node_pointers())
    {
      ++_errors;
      return nullptr;
    }
   
  if (_do_cache)
    {
      std::map<TrkrDefs::cluskey, PHG4Particle*>::iterator iter =
        _cache_max_truth_particle_by_energy.find(cluster_key);
      if (iter != _cache_max_truth_particle_by_energy.end())
	{
	  return iter->second;
	}
    }
  
  // loop over all particles associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(cluster_key);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
    {
      PHG4Particle* particle = *iter;
      float e = get_energy_contribution(cluster_key, particle);
      if (e > max_e)
	{
	  max_e = e;
	  max_particle = particle;
	}
    }
 
  if (_do_cache) _cache_max_truth_particle_by_energy.insert(make_pair(cluster_key, max_particle));
  
  return max_particle;
}

std::set<TrkrDefs::cluskey> SvtxClusterEval::all_clusters_from(PHG4Particle* truthparticle)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return std::set<TrkrDefs::cluskey>();
    }
  
  if (_strict)
    {
      assert(truthparticle);
    }
  else if (!truthparticle)
    {
      ++_errors;
      return std::set<TrkrDefs::cluskey>();
    }
  //check if cache is filled, if not fill it.
  //  if(_cache_all_clusters_from_particle.count(truthparticle)==0){
  if(_cache_all_clusters_from_particle.empty()){
    FillRecoClusterFromG4HitCache();
  }

  if (_do_cache)
    {
      std::map<PHG4Particle*, std::set<TrkrDefs::cluskey> >::iterator iter =
        _cache_all_clusters_from_particle.find(truthparticle);
      if (iter != _cache_all_clusters_from_particle.end())
	{
	  return iter->second;
	}
    }
  std::set<TrkrDefs::cluskey> clusters;
  return clusters;
}

void SvtxClusterEval::FillRecoClusterFromG4HitCache(){
  auto Mytimer = std::make_unique<PHTimer>("ReCl_timer");
  Mytimer->stop();
  Mytimer->restart();

  std::multimap<PHG4Particle*, TrkrDefs::cluskey> temp_clusters_from_particles;
  // loop over all the clusters
  for(const auto& hitsetkey:_clustermap->getHitSetKeys())
  {
    auto range = _clustermap->getClusters(hitsetkey);
    for( auto iter = range.first; iter != range.second; ++iter ){
      TrkrDefs::cluskey cluster_key = iter->first;
      
      // loop over all truth particles connected to this cluster
      std::set<PHG4Particle*> particles = all_truth_particles(cluster_key);
      for (std::set<PHG4Particle*>::iterator jter = particles.begin();
	   jter != particles.end();
	   ++jter){
	PHG4Particle* candidate = *jter;
	temp_clusters_from_particles.insert(make_pair(candidate, cluster_key));
      }
    }
  }
  //Loop over particles and fill cache
  PHG4TruthInfoContainer::ConstRange range = _truthinfo->GetParticleRange();
  for(PHG4TruthInfoContainer::ConstIterator iter = range.first;
      iter != range.second; ++iter){
    PHG4Particle* g4particle = iter->second;
    std::set<TrkrDefs::cluskey> clusters;
    std::multimap<PHG4Particle*, TrkrDefs::cluskey>::const_iterator lower_bound = temp_clusters_from_particles.lower_bound(g4particle);
    std::multimap<PHG4Particle*, TrkrDefs::cluskey>::const_iterator upper_bound = temp_clusters_from_particles.upper_bound(g4particle);
    std::multimap<PHG4Particle*, TrkrDefs::cluskey>::const_iterator cfp_iter;
    for(cfp_iter = lower_bound;cfp_iter != upper_bound;++cfp_iter){
      TrkrDefs::cluskey cluster_key = cfp_iter->second;
      clusters.insert(cluster_key);
    }
    _cache_all_clusters_from_particle.insert(make_pair(g4particle, clusters));
  }

  Mytimer->stop();

}

std::set<TrkrDefs::cluskey> SvtxClusterEval::all_clusters_from(PHG4Hit* truthhit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<TrkrDefs::cluskey>();
  }

  if (_strict)
  {
    assert(truthhit);
  }
  else if (!truthhit)
  {
    ++_errors;
    return std::set<TrkrDefs::cluskey>();
  }

  // one time, fill cache of g4hit/cluster pairs
  if(_cache_all_clusters_from_g4hit.size() == 0)
    {
      // make a map of truthhit, cluster_key inside this loop
      std::multimap<PHG4HitDefs::keytype, TrkrDefs::cluskey> truth_cluster_map;
      std::set<PHG4HitDefs::keytype> all_g4hits_set;
      std::map<PHG4HitDefs::keytype, PHG4Hit*> all_g4hits_map;

      // get all reco clusters
      if(_verbosity > 1) cout << "all_clusters_from_g4hit: list all reco clusters " << endl;
    
  for(const auto& hitsetkey:_clustermap->getHitSetKeys())
  {
    auto range = _clustermap->getClusters(hitsetkey);
	for( auto iter = range.first; iter != range.second; ++iter ){
	  TrkrDefs::cluskey cluster_key = iter->first;
	  int layer = TrkrDefs::getLayer(cluster_key);
	  TrkrCluster *clus = iter->second;
	  if(_verbosity > 1) 
	    {
	      cout << " layer " << layer << " cluster_key " << cluster_key << " adc " << clus->getAdc() 
		   << " localx " << clus->getLocalX() 
		   << " localy " << clus->getLocalY() 
		   << endl;
	      cout << "  associated hits:";
	      std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator> 
		hitrange = _cluster_hit_map->getHits(cluster_key);  // returns range of pairs {cluster key, hit key} for this cluskey
	      for(std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator
		    clushititer = hitrange.first; clushititer != hitrange.second; ++clushititer)
		{
		  TrkrDefs::hitkey hitkey = clushititer->second;
		  cout << " " << hitkey;
		}
	      cout << endl;
	    }
	  
	  // the returned truth hits were obtained from TrkrAssoc maps
	  std::set<PHG4Hit*> hits = all_truth_hits(cluster_key);
	  for (std::set<PHG4Hit*>::iterator jter = hits.begin();
	       jter != hits.end();
	       ++jter)
	    {
	      PHG4Hit* candidate = *jter;
	      PHG4HitDefs::keytype g4hit_key = candidate->get_hit_id();
	      
	      if(_verbosity > 5) 
		{
		  int gtrackID = candidate->get_trkid();
		  cout << "   adding cluster with cluster_key " << cluster_key << " g4hit with g4hit_key " << g4hit_key
		       << " gtrackID " << gtrackID 
		       << endl;
		}
	      
	      all_g4hits_set.insert(g4hit_key);
	      all_g4hits_map.insert(std::make_pair(g4hit_key, candidate));
	      
	      truth_cluster_map.insert(std::make_pair(g4hit_key, cluster_key));
	    }	  
	}
      }

      // now fill the cache
      // loop over all entries in all_g4hits
      for(std::set<PHG4HitDefs::keytype>::iterator iter = all_g4hits_set.begin(); iter != all_g4hits_set.end(); ++iter)
	{
	  PHG4HitDefs::keytype g4hit_key = *iter;
	  if(_verbosity > 5) cout << " associations for g4hit_key " << g4hit_key << endl;

	  std::map<PHG4HitDefs::keytype, PHG4Hit*>::iterator it = all_g4hits_map.find(g4hit_key);
 	  PHG4Hit *g4hit = it->second;

	  std::set<TrkrDefs::cluskey> assoc_clusters;

	  std::pair<std::multimap<PHG4HitDefs::keytype, TrkrDefs::cluskey>::iterator, 
		    std::multimap<PHG4HitDefs::keytype, TrkrDefs::cluskey>::iterator>  ret = truth_cluster_map.equal_range(g4hit_key);
	  for(std::multimap<PHG4HitDefs::keytype, TrkrDefs::cluskey>::iterator jter = ret.first; jter != ret.second; ++jter)
	    {
	      assoc_clusters.insert(jter->second);

	      if(_verbosity > 5)  cout << "             g4hit_key " << g4hit_key << " associated with cluster_key " << jter->second << endl; 
	    }
	  // done with this g4hit
	  _cache_all_clusters_from_g4hit.insert(make_pair(g4hit, assoc_clusters));
	}

    }
  
  // get the clusters
  std::set<TrkrDefs::cluskey> clusters;
  std::map<PHG4Hit*, std::set<TrkrDefs::cluskey> >::iterator iter =
    _cache_all_clusters_from_g4hit.find(truthhit);
  if (iter != _cache_all_clusters_from_g4hit.end())
    {
      return iter->second;
    }


  if (_clusters_per_layer.size() == 0)
  {
    fill_cluster_layer_map();
  }
      
   return clusters;
}

TrkrDefs::cluskey SvtxClusterEval::best_cluster_by_nhit(int gid, int layer)
{
   TrkrDefs::cluskey val =0;
  if (!has_node_pointers())
  {
    ++_errors;
    return  val;
  }

  /*  if (_strict)
  {
    assert(truthhit);
  }
  else if (!truthhit)
  {
    ++_errors;
    return val;
  }
  */
  // one time, fill cache of g4hit/cluster pairs
  if(_cache_best_cluster_from_gtrackid_layer.size() == 0){
    // get all reco clusters
    // cout << "cache size ==0" << endl;
    if(_verbosity > 1) cout << "all_clusters: found # " << _clustermap->size() << endl;
    
  for(const auto& hitsetkey:_clustermap->getHitSetKeys())
  {
      auto range = _clustermap->getClusters(hitsetkey);
      for( auto iter = range.first; iter != range.second; ++iter ){
	TrkrDefs::cluskey cluster_key = iter->first;
	int layer_in = TrkrDefs::getLayer(cluster_key);

	if(layer_in<0) continue;
	// TrkrCluster *clus = iter->second;
	
	std::pair<int, int> gid_lay = gtrackid_and_layer_by_nhit(cluster_key);

	//      std::map<std::pair<int, unsigned int>, TrkrDefs::cluskey>::iterator it_exists;
	//      it_exists = 
	if(_cache_best_cluster_from_gtrackid_layer.count(gid_lay)==0){
	  if(gid_lay.second >=0)
	    _cache_best_cluster_from_gtrackid_layer.insert(make_pair(gid_lay, cluster_key));
	}
	else
	  if(_verbosity > 2){ cout <<  "found doublematch" << endl;
	    cout << "ckey: " << cluster_key << " gtrackID: " << gid_lay.first << " layer: " << gid_lay.second << endl; 
	  }
      }
    }
  }
  

  // get the clusters
  TrkrDefs::cluskey best_cluster = 0;
  //  PHG4Hit*, std::set<TrkrDefs::cluskey> >::iterator iter =

  std::map<std::pair<int, int>, TrkrDefs::cluskey>::iterator iter = _cache_best_cluster_from_gtrackid_layer.find(make_pair(gid,layer));
  if (iter != _cache_best_cluster_from_gtrackid_layer.end())
    {
      return iter->second;
    }

  return best_cluster;
}

TrkrDefs::cluskey SvtxClusterEval::best_cluster_from(PHG4Hit* truthhit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return 0;
  }

  if (_strict)
  {
    assert(truthhit);
  }
  else if (!truthhit)
  {
    ++_errors;
    return 0;
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, TrkrDefs::cluskey>::iterator iter =
        _cache_best_cluster_from_g4hit.find(truthhit);
    if (iter != _cache_best_cluster_from_g4hit.end())
    {
      return iter->second;
    }
  }

  TrkrDefs::cluskey best_cluster = 0;
  float best_energy = 0.0;
  std::set<TrkrDefs::cluskey> clusters = all_clusters_from(truthhit);
  for (std::set<TrkrDefs::cluskey>::iterator iter = clusters.begin();
       iter != clusters.end();
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = *iter;
    float energy = get_energy_contribution(cluster_key, truthhit);
    if (energy > best_energy)
    {
      best_cluster = cluster_key;
      best_energy = energy;
    }
  }

  if (_do_cache) _cache_best_cluster_from_g4hit.insert(make_pair(truthhit, best_cluster));

  return best_cluster;
}

// overlap calculations
float SvtxClusterEval::get_energy_contribution(TrkrDefs::cluskey cluster_key, PHG4Particle* particle)
{
  // Note: this does not work correctly for the TPC
  // It assumes one g4hit per layer. Use the truth cluster energy instead.

  if (!has_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
//    assert(cluster_key);
    assert(particle);
  }
  else if ( !particle)
  {
    ++_errors;
    return NAN;
  }

  if (_do_cache)
  {
    std::map<std::pair<TrkrDefs::cluskey, PHG4Particle*>, float>::iterator iter =
        _cache_get_energy_contribution_g4particle.find(make_pair(cluster_key, particle));
    if (iter != _cache_get_energy_contribution_g4particle.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;
  std::set<PHG4Hit*> hits = all_truth_hits(cluster_key);
  for (std::set<PHG4Hit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter)
  {
    PHG4Hit* hit = *iter;
    if (get_truth_eval()->is_g4hit_from_particle(hit, particle))
    {
      energy += hit->get_edep();
    }
  }

  if (_do_cache) _cache_get_energy_contribution_g4particle.insert(make_pair(make_pair(cluster_key, particle), energy));

  return energy;
}

float SvtxClusterEval::get_energy_contribution(TrkrDefs::cluskey cluster_key, PHG4Hit* g4hit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
//    assert(cluster_key);
    assert(g4hit);
  }
  else if ( !g4hit)
  {
    ++_errors;
    return NAN;
  }

  if ((_do_cache) &&
      (_cache_get_energy_contribution_g4hit.find(make_pair(cluster_key, g4hit)) !=
       _cache_get_energy_contribution_g4hit.end()))
  {
    return _cache_get_energy_contribution_g4hit[make_pair(cluster_key, g4hit)];
  }

  // this is a fairly simple existance check right now, but might be more
  // complex in the future, so this is here mostly as future-proofing.

  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster_key);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter)
  {
    PHG4Hit* candidate = *iter;
    if (candidate->get_hit_id() != g4hit->get_hit_id()) continue;
    energy += candidate->get_edep();
  }

  if (_do_cache) _cache_get_energy_contribution_g4hit.insert(make_pair(make_pair(cluster_key, g4hit), energy));

  return energy;
}

void SvtxClusterEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...

  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if(!_clustermap) {
    _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  } 
  else 
    {
      /// Need a catch for if the node corrected cluster node exists but hasn't been filled 
      /// yet, in e.g. the case of truth track seeding
      if(_clustermap->size() == 0)
	_clustermap = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
    }
  

  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  _hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  _g4hits_mms = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");
  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  
  return;
}

void SvtxClusterEval::fill_cluster_layer_map()
{
  // loop over all the clusters
  for(const auto& hitsetkey:_clustermap->getHitSetKeys())
  {
    auto range = _clustermap->getClusters(hitsetkey);
    for( auto iter = range.first; iter != range.second; ++iter ){
      TrkrDefs::cluskey cluster_key = iter->first;
      unsigned int ilayer = TrkrDefs::getLayer(cluster_key);
      TrkrCluster *cluster = iter->second;
      auto glob = getGlobalPosition(cluster_key, cluster);
      float clus_phi = atan2(glob(1), glob(0));
      
      multimap<unsigned int, innerMap>::iterator it = _clusters_per_layer.find(ilayer);
      if (it == _clusters_per_layer.end())
	{
	  it = _clusters_per_layer.insert(make_pair(ilayer, innerMap()));
	}
      it->second.insert(make_pair(clus_phi, cluster_key));
      
      //address wrapping along +/-PI by filling larger area of the map
      if (clus_phi - (-M_PI) < _clusters_searching_window) it->second.insert(make_pair(clus_phi + 2 * M_PI, cluster_key));
      if (M_PI - clus_phi < _clusters_searching_window) it->second.insert(make_pair(clus_phi - 2 * M_PI, cluster_key));
    }
  }
  return;
}

bool SvtxClusterEval::has_node_pointers()
{
  if (_strict)
    assert(_clustermap);
  else if (!_clustermap)
    return false;

  if (_strict)
    assert(_truthinfo);
  else if (!_truthinfo)
    return false;

  return true;
}

float SvtxClusterEval::fast_approx_atan2(float y, float x)
{
  if (x != 0.0f)
  {
    if (fabsf(x) > fabsf(y))
    {
      const float z = y / x;
      if (x > 0.0)
      {
        // atan2(y,x) = atan(y/x) if x > 0
        return fast_approx_atan2(z);
      }
      else if (y >= 0.0)
      {
        // atan2(y,x) = atan(y/x) + PI if x < 0, y >= 0
        return fast_approx_atan2(z) + M_PI;
      }
      else
      {
        // atan2(y,x) = atan(y/x) - PI if x < 0, y < 0
        return fast_approx_atan2(z) - M_PI;
      }
    }
    else  // Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
    {
      const float z = x / y;
      if (y > 0.0)
      {
        // atan2(y,x) = PI/2 - atan(x/y) if |y/x| > 1, y > 0
        return -fast_approx_atan2(z) + M_PI_2;
      }
      else
      {
        // atan2(y,x) = -PI/2 - atan(x/y) if |y/x| > 1, y < 0
        return -fast_approx_atan2(z) - M_PI_2;
      }
    }
  }
  else
  {
    if (y > 0.0f)  // x = 0, y > 0
    {
      return M_PI_2;
    }
    else if (y < 0.0f)  // x = 0, y < 0
    {
      return -M_PI_2;
    }
  }
  return 0.0f;  // x,y = 0. Could return NaN instead.
}

float SvtxClusterEval::fast_approx_atan2(float z)
{
  // Polynomial approximating arctangenet on the range -1,1.
  // Max error < 0.005 (or 0.29 degrees)
  const float n1 = 0.97239411f;
  const float n2 = -0.19194795f;
  return (n1 + n2 * z * z) * z;
}

Acts::Vector3 SvtxClusterEval::getGlobalPosition(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster)
{
  Acts::Vector3 glob;
  glob = _tgeometry->getGlobalPosition(cluster_key, cluster);

  return glob;
}
