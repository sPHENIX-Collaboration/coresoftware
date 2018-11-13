#include "SvtxClusterEval.h"

#include "SvtxHitEval.h"

#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxHit.h>
#include <g4hough/SvtxHitMap.h>
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
  , _hitmap(nullptr)
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

std::set<PHG4Hit*> SvtxClusterEval::all_truth_hits(SvtxCluster* cluster)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return std::set<PHG4Hit*>();
  }

  if (_do_cache)
  {
    std::map<SvtxCluster*, std::set<PHG4Hit*> >::iterator iter =
        _cache_all_truth_hits.find(cluster);
    if (iter != _cache_all_truth_hits.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all hit cells
  for (SvtxCluster::ConstHitIter hiter = cluster->begin_hits();
       hiter != cluster->end_hits();
       ++hiter)
  {
    SvtxHit* hit = _hitmap->get(*hiter);

    if (_strict)
    {
      assert(hit);
    }
    else if (!hit)
    {
      ++_errors;
      continue;
    }

    std::set<PHG4Hit*> new_g4hits = _hiteval.all_truth_hits(hit);

    for (std::set<PHG4Hit*>::iterator iter = new_g4hits.begin();
         iter != new_g4hits.end();
         ++iter)
    {
      //cout << "cluster_eval insert g4hit " << *iter->
      truth_hits.insert(*iter);
    }
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(cluster, truth_hits));

  return truth_hits;
}

PHG4Hit* SvtxClusterEval::max_truth_hit_by_energy(SvtxCluster* cluster)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<SvtxCluster*, PHG4Hit*>::iterator iter =
        _cache_max_truth_hit_by_energy.find(cluster);
    if (iter != _cache_max_truth_hit_by_energy.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> hits = all_truth_hits(cluster);
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

  if (_do_cache) _cache_max_truth_hit_by_energy.insert(make_pair(cluster, max_hit));

  return max_hit;
}

std::set<PHG4Particle*> SvtxClusterEval::all_truth_particles(SvtxCluster* cluster)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return std::set<PHG4Particle*>();
  }

  if (_do_cache)
  {
    std::map<SvtxCluster*, std::set<PHG4Particle*> >::iterator iter =
        _cache_all_truth_particles.find(cluster);
    if (iter != _cache_all_truth_particles.end())
    {
      return iter->second;
    }
  }

  std::set<PHG4Particle*> truth_particles;

  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster);

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

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(cluster, truth_particles));

  return truth_particles;
}

PHG4Particle* SvtxClusterEval::max_truth_particle_by_energy(SvtxCluster* cluster)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(cluster);
  }
  else if (!cluster)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<SvtxCluster*, PHG4Particle*>::iterator iter =
        _cache_max_truth_particle_by_energy.find(cluster);
    if (iter != _cache_max_truth_particle_by_energy.end())
    {
      return iter->second;
    }
  }

  // loop over all particles associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = nullptr;
  float max_e = FLT_MAX * -1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(cluster);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter)
  {
    PHG4Particle* particle = *iter;
    float e = get_energy_contribution(cluster, particle);
    if (e > max_e)
    {
      max_e = e;
      max_particle = particle;
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_energy.insert(make_pair(cluster, max_particle));

  return max_particle;
}

std::set<SvtxCluster*> SvtxClusterEval::all_clusters_from(PHG4Particle* truthparticle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<SvtxCluster*>();
  }

  if (_strict)
  {
    assert(truthparticle);
  }
  else if (!truthparticle)
  {
    ++_errors;
    return std::set<SvtxCluster*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Particle*, std::set<SvtxCluster*> >::iterator iter =
        _cache_all_clusters_from_particle.find(truthparticle);
    if (iter != _cache_all_clusters_from_particle.end())
    {
      return iter->second;
    }
  }

  std::set<SvtxCluster*> clusters;

  // loop over all the clusters
  for (SvtxClusterMap::Iter iter = _clustermap->begin();
       iter != _clustermap->end();
       ++iter)
  {
    SvtxCluster* cluster = iter->second;

    // loop over all truth particles connected to this cluster
    std::set<PHG4Particle*> particles = all_truth_particles(cluster);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
         jter != particles.end();
         ++jter)
    {
      PHG4Particle* candidate = *jter;
      if (get_truth_eval()->are_same_particle(candidate, truthparticle))
      {
        clusters.insert(cluster);
      }
    }
  }

  if (_do_cache) _cache_all_clusters_from_particle.insert(make_pair(truthparticle, clusters));

  return clusters;
}

std::set<SvtxCluster*> SvtxClusterEval::all_clusters_from(PHG4Hit* truthhit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return std::set<SvtxCluster*>();
  }

  if (_strict)
  {
    assert(truthhit);
  }
  else if (!truthhit)
  {
    ++_errors;
    return std::set<SvtxCluster*>();
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, std::set<SvtxCluster*> >::iterator iter =
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
  std::set<SvtxCluster*> clusters;

  unsigned int hit_layer = truthhit->get_layer();
  // loop over all the clusters

  int count = 0;
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

    for (multimap<float, SvtxCluster*>::iterator liter = iter_lower_bound;
         liter != iter_upper_bound;
         liter++)
    {
      count++;
      /*
      for (SvtxClusterMap::Iter iter = _clustermap->begin();
      iter != _clustermap->end();
      ++iter) {
      
      //
      SvtxCluster* cluster = iter->second;
    */
      SvtxCluster* cluster = liter->second;
      if (cluster->get_layer() != hit_layer) continue;

      // loop over all truth hits connected to this cluster
      std::set<PHG4Hit*> hits = all_truth_hits(cluster);
      for (std::set<PHG4Hit*>::iterator jter = hits.begin();
           jter != hits.end();
           ++jter)
      {
        PHG4Hit* candidate = *jter;
        if (candidate->get_hit_id() == truthhit->get_hit_id())
        {
          clusters.insert(cluster);
        }
      }  //      for (std::set<PHG4Hit*>::iterator jter = hits.begin();

    }  //    for (multimap<float, SvtxCluster*>::iterator liter = iter_lower_bound;

  }  //  if (miter != _clusters_per_layer.end())

  //  cout << "count " << count << endl;
  if (_do_cache) _cache_all_clusters_from_g4hit.insert(make_pair(truthhit, clusters));

  return clusters;
}

SvtxCluster* SvtxClusterEval::best_cluster_from(PHG4Hit* truthhit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return nullptr;
  }

  if (_strict)
  {
    assert(truthhit);
  }
  else if (!truthhit)
  {
    ++_errors;
    return nullptr;
  }

  if (_do_cache)
  {
    std::map<PHG4Hit*, SvtxCluster*>::iterator iter =
        _cache_best_cluster_from_g4hit.find(truthhit);
    if (iter != _cache_best_cluster_from_g4hit.end())
    {
      return iter->second;
    }
  }

  SvtxCluster* best_cluster = nullptr;
  float best_energy = 0.0;
  std::set<SvtxCluster*> clusters = all_clusters_from(truthhit);
  for (std::set<SvtxCluster*>::iterator iter = clusters.begin();
       iter != clusters.end();
       ++iter)
  {
    SvtxCluster* cluster = *iter;
    float energy = get_energy_contribution(cluster, truthhit);
    if (energy > best_energy)
    {
      best_cluster = cluster;
      best_energy = energy;
    }
  }

  if (_do_cache) _cache_best_cluster_from_g4hit.insert(make_pair(truthhit, best_cluster));

  return best_cluster;
}

// overlap calculations
float SvtxClusterEval::get_energy_contribution(SvtxCluster* cluster, PHG4Particle* particle)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
    assert(cluster);
    assert(particle);
  }
  else if (!cluster || !particle)
  {
    ++_errors;
    return NAN;
  }

  if (_do_cache)
  {
    std::map<std::pair<SvtxCluster*, PHG4Particle*>, float>::iterator iter =
        _cache_get_energy_contribution_g4particle.find(make_pair(cluster, particle));
    if (iter != _cache_get_energy_contribution_g4particle.end())
    {
      return iter->second;
    }
  }

  float energy = 0.0;
  std::set<PHG4Hit*> hits = all_truth_hits(cluster);
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

  if (_do_cache) _cache_get_energy_contribution_g4particle.insert(make_pair(make_pair(cluster, particle), energy));

  return energy;
}

float SvtxClusterEval::get_energy_contribution(SvtxCluster* cluster, PHG4Hit* g4hit)
{
  if (!has_node_pointers())
  {
    ++_errors;
    return NAN;
  }

  if (_strict)
  {
    assert(cluster);
    assert(g4hit);
  }
  else if (!cluster || !g4hit)
  {
    ++_errors;
    return NAN;
  }

  if ((_do_cache) &&
      (_cache_get_energy_contribution_g4hit.find(make_pair(cluster, g4hit)) !=
       _cache_get_energy_contribution_g4hit.end()))
  {
    return _cache_get_energy_contribution_g4hit[make_pair(cluster, g4hit)];
  }

  // this is a fairly simple existance check right now, but might be more
  // complex in the future, so this is here mostly as future-proofing.

  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter)
  {
    PHG4Hit* candidate = *iter;
    if (candidate->get_hit_id() != g4hit->get_hit_id()) continue;
    energy += candidate->get_edep();
  }

  if (_do_cache) _cache_get_energy_contribution_g4hit.insert(make_pair(make_pair(cluster, g4hit), energy));

  return energy;
}

void SvtxClusterEval::get_node_pointers(PHCompositeNode* topNode)
{
  // need things off of the DST...
  _clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");

  _hitmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  return;
}

void SvtxClusterEval::fill_cluster_layer_map()
{
  // loop over all the clusters

  //  for(unsigned int i = 0;i<47;i++){
  //    _clusters_per_layer.insert(make_pair(i,innerMap()));
  //  }

  for (SvtxClusterMap::Iter iter = _clustermap->begin(); iter != _clustermap->end(); ++iter)
  {
    SvtxCluster* cluster = iter->second;
    unsigned int ilayer = cluster->get_layer();
    float clus_phi = fast_approx_atan2(cluster->get_y(), cluster->get_x());

    multimap<unsigned int, innerMap>::iterator it = _clusters_per_layer.find(ilayer);
    if (it == _clusters_per_layer.end())
    {
      it = _clusters_per_layer.insert(make_pair(ilayer, innerMap()));
    }
    it->second.insert(make_pair(clus_phi, cluster));

    //address wrapping along +/-PI by filling larger area of the map
    if (clus_phi - (-M_PI) < _clusters_searching_window) it->second.insert(make_pair(clus_phi + 2 * M_PI, cluster));
    if (M_PI - clus_phi < _clusters_searching_window) it->second.insert(make_pair(clus_phi - 2 * M_PI, cluster));
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
    assert(_hitmap);
  else if (!_hitmap)
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
