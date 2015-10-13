
#include "SvtxClusterEval.h"

#include "SvtxHitEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <cstdlib>
#include <set>
#include <map>
#include <float.h>
#include <algorithm>
#include <cassert>

using namespace std;

SvtxClusterEval::SvtxClusterEval(PHCompositeNode* topNode)
  : _hiteval(topNode),
    _clustermap(NULL),
    _hitmap(NULL),
    _truthinfo(NULL),
    _strict(true),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_max_truth_hit_by_energy(),
    _cache_all_truth_particles(),
    _cache_max_truth_particle_by_energy(),
    _cache_all_clusters_from_particle(),
    _cache_all_clusters_from_g4hit(),
    _cache_best_cluster_from_g4hit(),
    _cache_get_energy_contribution_g4particle(),
    _cache_get_energy_contribution_g4hit() {
  get_node_pointers(topNode);
}

void SvtxClusterEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_max_truth_hit_by_energy.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_energy.clear();
  _cache_all_clusters_from_particle.clear();
  _cache_all_clusters_from_g4hit.clear();
  _cache_best_cluster_from_g4hit.clear();
  _cache_get_energy_contribution_g4particle.clear();
  _cache_get_energy_contribution_g4hit.clear();

  _hiteval.next_event(topNode);
  
  get_node_pointers(topNode);
}

std::set<PHG4Hit*> SvtxClusterEval::all_truth_hits(SvtxCluster* cluster) {

  if (_strict) assert(cluster);
  else if (!cluster) return std::set<PHG4Hit*>();
  
  if (_do_cache) {
    std::map<SvtxCluster*,std::set<PHG4Hit*> >::iterator iter =
      _cache_all_truth_hits.find(cluster);
    if (iter != _cache_all_truth_hits.end()) {
      return iter->second;
    }
  }
 
  std::set<PHG4Hit*> truth_hits;
  
  // loop over all hit cells
  for (SvtxCluster::ConstHitIter hiter = cluster->begin_hits();
       hiter != cluster->end_hits();
       ++hiter) {
    SvtxHit* hit = _hitmap->get(*hiter);

    if (_strict) assert(hit);
    else if (!hit) continue;
    
    std::set<PHG4Hit*> new_g4hits = _hiteval.all_truth_hits(hit);

    for (std::set<PHG4Hit*>::iterator iter = new_g4hits.begin();
	 iter != new_g4hits.end();
	 ++iter) {

      truth_hits.insert(*iter);
    }
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(cluster,truth_hits));
  
  return truth_hits;
}

PHG4Hit* SvtxClusterEval::max_truth_hit_by_energy(SvtxCluster* cluster) {

  if (_strict) assert(cluster);
  else if (!cluster) return NULL;
  
  if (_do_cache) {
    std::map<SvtxCluster*,PHG4Hit*>::iterator iter =
      _cache_max_truth_hit_by_energy.find(cluster);
    if (iter != _cache_max_truth_hit_by_energy.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> hits = all_truth_hits(cluster);
  PHG4Hit* max_hit = NULL;
  float max_e = FLT_MAX*-1.0;
  for (std::set<PHG4Hit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter) {
    PHG4Hit *hit = *iter;
    if (hit->get_edep() > max_e) {
      max_e = hit->get_edep();
      max_hit = hit;
    }
  }

  if (_do_cache) _cache_max_truth_hit_by_energy.insert(make_pair(cluster,max_hit));
  
  return max_hit;
}
  
std::set<PHG4Particle*> SvtxClusterEval::all_truth_particles(SvtxCluster* cluster) {

  if (_strict) assert(cluster);
  else if (!cluster) return std::set<PHG4Particle*>();
  
  if (_do_cache) {
    std::map<SvtxCluster*,std::set<PHG4Particle*> >::iterator iter =
      _cache_all_truth_particles.find(cluster);
    if (iter !=	_cache_all_truth_particles.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Particle*> truth_particles;
  
  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster);

  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* hit = *iter;
    PHG4Particle* particle = _truthinfo->GetHit( hit->get_trkid() );

    if (_strict) assert(particle);
    else if (!particle) continue;
    
    truth_particles.insert(particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(cluster,truth_particles));
  
  return truth_particles;
}

PHG4Particle* SvtxClusterEval::max_truth_particle_by_energy(SvtxCluster* cluster) {

  if (_strict) assert(cluster);
  else if (!cluster) return NULL;
  
  if (_do_cache) {
    std::map<SvtxCluster*,PHG4Particle*>::iterator iter =
      _cache_max_truth_particle_by_energy.find(cluster);
    if (iter != _cache_max_truth_particle_by_energy.end()) {
      return iter->second;
    }
  }

  // loop over all particles associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = NULL;
  float max_e = FLT_MAX*-1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(cluster);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter) {

    PHG4Particle* particle = *iter;
    float e = get_energy_contribution(cluster,particle);
    if (e > max_e) {
      max_e = e;
      max_particle = particle;      
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_energy.insert(make_pair(cluster,max_particle));
  
  return max_particle;
}

std::set<SvtxCluster*> SvtxClusterEval::all_clusters_from(PHG4Particle* truthparticle) { 

  if (_strict) assert(truthparticle);
  else if (!truthparticle) return std::set<SvtxCluster*>();
  
  if (_do_cache) {
    std::map<PHG4Particle*,std::set<SvtxCluster*> >::iterator iter =
      _cache_all_clusters_from_particle.find(truthparticle);
    if (iter != _cache_all_clusters_from_particle.end()) {
      return iter->second;
    }
  }
  
  std::set<SvtxCluster*> clusters;
  
  // loop over all the clusters
  for (SvtxClusterMap::Iter iter = _clustermap->begin();
       iter != _clustermap->end();
       ++iter) {

    SvtxCluster* cluster = iter->second;

    // loop over all truth particles connected to this cluster
    std::set<PHG4Particle*> particles = all_truth_particles(cluster);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
	 jter != particles.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;
      if (candidate->get_track_id() == truthparticle->get_track_id()) {
	clusters.insert(cluster);
      }    
    }
  }

  if (_do_cache) _cache_all_clusters_from_particle.insert(make_pair(truthparticle,clusters));
  
  return clusters;
}

std::set<SvtxCluster*> SvtxClusterEval::all_clusters_from(PHG4Hit* truthhit) {

  if (_strict) assert(truthhit);
  else if (!truthhit) return std::set<SvtxCluster*>();
  
  if (_do_cache) {
    std::map<PHG4Hit*,std::set<SvtxCluster*> >::iterator iter =
      _cache_all_clusters_from_g4hit.find(truthhit);
    if (iter != _cache_all_clusters_from_g4hit.end()) {
      return iter->second;
    }
  }
  
  std::set<SvtxCluster*> clusters;
  
  // loop over all the clusters
  for (SvtxClusterMap::Iter iter = _clustermap->begin();
       iter != _clustermap->end();
       ++iter) {

    SvtxCluster* cluster = iter->second;

    // loop over all truth hits connected to this cluster
    std::set<PHG4Hit*> hits = all_truth_hits(cluster);
    for (std::set<PHG4Hit*>::iterator jter = hits.begin();
	 jter != hits.end();
	 ++jter) {
      PHG4Hit* candidate = *jter;
      if (candidate->get_hit_id() == truthhit->get_hit_id()) {
	clusters.insert(cluster);
      }    
    }
  }

  if (_do_cache) _cache_all_clusters_from_g4hit.insert(make_pair(truthhit,clusters));
  
  return clusters;
}

SvtxCluster* SvtxClusterEval::best_cluster_from(PHG4Hit* truthhit) {

  if (_strict) assert(truthhit);
  else if (!truthhit) return NULL;
  
  if (_do_cache) {
    std::map<PHG4Hit*,SvtxCluster*>::iterator iter =
      _cache_best_cluster_from_g4hit.find(truthhit);
    if (iter != _cache_best_cluster_from_g4hit.end()) {
      return iter->second;
    }
  }

  SvtxCluster* best_cluster = NULL;
  float best_energy = 0.0;  
  std::set<SvtxCluster*> clusters = all_clusters_from(truthhit);
  for (std::set<SvtxCluster*>::iterator iter = clusters.begin();
       iter != clusters.end();
       ++iter) {
    SvtxCluster* cluster = *iter;
    float energy = get_energy_contribution(cluster,truthhit);
    if (energy > best_energy) {
      best_cluster = cluster;
      best_energy = energy;
    }
  }
 
  if (_do_cache) _cache_best_cluster_from_g4hit.insert(make_pair(truthhit,best_cluster));
  
  return best_cluster;
}
  
// overlap calculations
float SvtxClusterEval::get_energy_contribution(SvtxCluster* cluster, PHG4Particle* particle) {

  if (_strict) {
    assert(cluster);
    assert(particle);
  } else if (!cluster||!particle) {
    return NAN;
  }

  if (_do_cache) {
    std::map<std::pair<SvtxCluster*,PHG4Particle*>, float>::iterator iter =
      _cache_get_energy_contribution_g4particle.find(make_pair(cluster,particle));
    if (iter !=	_cache_get_energy_contribution_g4particle.end()) {
      return iter->second;
    }
  }
  
  float energy = 0.0;
  std::set<PHG4Hit*> hits = all_truth_hits(cluster);
  for (std::set<PHG4Hit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter) {
    PHG4Hit* hit = *iter;
    if (hit->get_trkid() == particle->get_track_id()) {
      energy += hit->get_edep();
    }
  }

  if (_do_cache) _cache_get_energy_contribution_g4particle.insert(make_pair(make_pair(cluster,particle),energy));
  
  return energy;
}

float SvtxClusterEval::get_energy_contribution(SvtxCluster* cluster, PHG4Hit* g4hit) {

 if (_strict) {
    assert(cluster);
    assert(g4hit);
  } else if (!cluster||!g4hit) {
    return NAN;
  }
  
  if ((_do_cache) &&
      (_cache_get_energy_contribution_g4hit.find(make_pair(cluster,g4hit)) !=
       _cache_get_energy_contribution_g4hit.end())) {
    return _cache_get_energy_contribution_g4hit[make_pair(cluster,g4hit)];
  }

  // this is a fairly simple existance check right now, but might be more
  // complex in the future, so this is here mostly as future-proofing.
  
  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* candidate = *iter;
    if (candidate->get_hit_id() != g4hit->get_hit_id()) continue;    
    energy += candidate->get_edep();
  }

  if (_do_cache) _cache_get_energy_contribution_g4hit.insert(make_pair(make_pair(cluster,g4hit),energy));
  
  return energy;
}

void SvtxClusterEval::get_node_pointers(PHCompositeNode *topNode) {

  // need things off of the DST...
  _clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if (!_clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }

  _hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hitmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxHitMap" << endl;
    exit(-1);
  }

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }
  
  return;
}
