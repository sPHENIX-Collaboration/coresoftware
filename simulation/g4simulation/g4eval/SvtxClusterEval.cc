#include "SvtxClusterEval.h"

#include "SvtxHitEval.h"

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TMath.h>

#include <TVector3.h>
#include <float.h>
#include <cassert>
#include <cmath>
#include <map>
#include <set>

using namespace std;

SvtxClusterEval::SvtxClusterEval(PHCompositeNode* topNode)
  : _hiteval(topNode)
  , _clustermap(nullptr)
    //  , _hitmap(nullptr)
  , _truthinfo(nullptr)
  , _strict(false)
  , _verbosity(1)
  , _errors(0)
  , _do_cache(true)
  , _cache_all_truth_hits()
  , _cache_max_truth_hit_by_energy()
  , _cache_all_truth_particles()
  , _cache_max_truth_particle_by_energy()
  , _cache_all_clusters_from_particle()
  , _cache_all_clusters_from_g4hit()
  , _cache_best_cluster_from_g4hit()
  , _cache_get_energy_contribution_g4particle()
  , _cache_get_energy_contribution_g4hit()
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
  _cache_max_truth_hit_by_energy.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_energy.clear();
  _cache_all_clusters_from_particle.clear();
  _cache_all_clusters_from_g4hit.clear();
  _cache_best_cluster_from_g4hit.clear();
  _cache_get_energy_contribution_g4particle.clear();
  _cache_get_energy_contribution_g4hit.clear();
  
  _clusters_per_layer.clear();
  //  _g4hits_per_layer.clear();
  _hiteval.next_event(topNode);
  
  get_node_pointers(topNode);
} 

std::set<PHG4Hit*> SvtxClusterEval::all_truth_hits(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return std::set<PHG4Hit*>();
    }
  
  if (_strict)
    {
      assert(cluster_key);
    }
  else if (!cluster_key)
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

	  // extract the g4 hit key here and add the hits to the set
	  PHG4HitDefs::keytype g4hitkey = htiter->second.second;
	  PHG4Hit * g4hit;
	  unsigned int trkrid = TrkrDefs::getTrkrId(hitsetkey);
	  if(trkrid == TrkrDefs::tpcId)
	    g4hit = _g4hits_tpc->findHit(g4hitkey);
	  else if(trkrid == TrkrDefs::inttId)
	    g4hit = _g4hits_intt->findHit(g4hitkey);
	  else
	    g4hit = _g4hits_mvtx->findHit(g4hitkey);
	  truth_hits.insert(g4hit);	      
	} // end loop over g4hits associated with hitsetkey and hitkey
    } // end loop over hits associated with cluskey  

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(cluster_key, truth_hits));

  return truth_hits;
}

PHG4Hit* SvtxClusterEval::max_truth_hit_by_energy(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return nullptr;
    }
  
  if (_strict)
    {
      assert(cluster_key);
    }
  else if (!cluster_key)
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
  
  if (_strict)
    {
      assert(cluster_key);
    }
  else if (!cluster_key)
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

PHG4Particle* SvtxClusterEval::max_truth_particle_by_energy(TrkrDefs::cluskey cluster_key)
{
  if (!has_node_pointers())
    {
      ++_errors;
      return nullptr;
    }
  
  if (_strict)
    {
      assert(cluster_key);
    }
  else if (!cluster_key)
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

  // loop over all the clusters
  TrkrClusterContainer::ConstRange all_clusters = _clustermap->getClusters();
  for (TrkrClusterContainer::ConstIterator iter = all_clusters.first;
       iter != all_clusters.second;
       ++iter)
  {
    TrkrDefs::cluskey cluster_key = iter->first;
    
    // loop over all truth particles connected to this cluster
    std::set<PHG4Particle*> particles = all_truth_particles(cluster_key);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
         jter != particles.end();
         ++jter)
    {
      PHG4Particle* candidate = *jter;
      if (get_truth_eval()->are_same_particle(candidate, truthparticle))
      {
        clusters.insert(cluster_key);
      }
    }
  }

  if (_do_cache) _cache_all_clusters_from_particle.insert(make_pair(truthparticle, clusters));

  return clusters;
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

  if (_do_cache)
  {
    std::map<PHG4Hit*, std::set<TrkrDefs::cluskey> >::iterator iter =
        _cache_all_clusters_from_g4hit.find(truthhit);
    if (iter != _cache_all_clusters_from_g4hit.end())
    {
      return iter->second;
    }
  }
  if (_clusters_per_layer.size() == 0)
  {
    fill_cluster_layer_map();
  }
  std::set<TrkrDefs::cluskey> clusters;

  unsigned int hit_layer = truthhit->get_layer();
  // loop over all the clusters

  multimap<unsigned int, innerMap>::iterator miter = _clusters_per_layer.find(hit_layer);
  if (miter != _clusters_per_layer.end())
  {
    const float hit_phi = fast_approx_atan2(truthhit->get_avg_y(), truthhit->get_avg_x());

    if (_verbosity >= 2)
    {
      cout << "SvtxClusterEval::all_clusters_from - hit_phi = " << hit_phi
           << ", miter->first = " << miter->first
           << ", clusters_searching_window = " << _clusters_searching_window
           << ", miter->second.size() = " << miter->second.size()
           << endl;
    }

    auto iter_lower_bound = miter->second.lower_bound(hit_phi - _clusters_searching_window);
    auto iter_upper_bound = miter->second.upper_bound(hit_phi + _clusters_searching_window);

    for (multimap<float, TrkrDefs::cluskey>::iterator liter = iter_lower_bound;
         liter != iter_upper_bound;
         ++liter)
    {
      TrkrDefs::cluskey cluster_key = liter->second;

      if (TrkrDefs::getLayer(cluster_key) != hit_layer) continue;

      // loop over all truth hits connected to this cluster
      std::set<PHG4Hit*> hits = all_truth_hits(cluster_key);
      for (std::set<PHG4Hit*>::iterator jter = hits.begin();
           jter != hits.end();
           ++jter)
      {
        PHG4Hit* candidate = *jter;
        if (candidate->get_hit_id() == truthhit->get_hit_id())
        {
          clusters.insert(cluster_key);
        }
      }  //      for (std::set<PHG4Hit*>::iterator jter = hits.begin();

    }  //    for (multimap<float, TrkrDefs::cluskey>::iterator liter = iter_lower_bound;

  }  //  if (miter != _clusters_per_layer.end())

  if (_do_cache) _cache_all_clusters_from_g4hit.insert(make_pair(truthhit, clusters));

  return clusters;
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
  if (!has_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
    assert(cluster_key);
    assert(particle);
  }
  else if (!cluster_key || !particle)
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
    assert(cluster_key);
    assert(g4hit);
  }
  else if (!cluster_key || !g4hit)
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

  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  _hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");

  return;
}

void SvtxClusterEval::fill_cluster_layer_map()
{
  // loop over all the clusters
  TrkrClusterContainer::ConstRange all_clusters = _clustermap->getClusters();
  for (TrkrClusterContainer::ConstIterator iter = all_clusters.first; iter != all_clusters.second; ++iter)
  {
    TrkrDefs::cluskey cluster_key = iter->first;
    unsigned int ilayer = TrkrDefs::getLayer(cluster_key);
    TrkrCluster *cluster = iter->second;
    float clus_phi = fast_approx_atan2(cluster->getY(), cluster->getX());

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
