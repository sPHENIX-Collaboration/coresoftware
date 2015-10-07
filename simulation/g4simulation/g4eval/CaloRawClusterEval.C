
#include "CaloRawClusterEval.h"
#include "CaloTruthEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <string>
#include <cstdlib>
#include <set>
#include <map>
#include <float.h>
#include <algorithm>
#include <cassert>

using namespace std;

CaloRawClusterEval::CaloRawClusterEval(PHCompositeNode* topNode, std::string caloname)
  : _caloname(caloname),
    _towereval(topNode,caloname),
    _clusters(NULL),
    _towers(NULL),
    _strict(true),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_primaries(),
    _cache_max_truth_primary_by_energy(),
    _cache_all_clusters_from_primary(),
    _cache_best_cluster_from_primary(),
    _cache_get_energy_contribution_primary() {
  get_node_pointers(topNode);
}

void CaloRawClusterEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_primaries.clear();
  _cache_max_truth_primary_by_energy.clear();
  _cache_all_clusters_from_primary.clear();
  _cache_best_cluster_from_primary.clear();
  _cache_get_energy_contribution_primary.clear();

  _towereval.next_event(topNode);
  
  get_node_pointers(topNode);
}

std::set<PHG4Hit*> CaloRawClusterEval::all_truth_hits(RawCluster* cluster) {

  if (_strict) assert(cluster);
  else if (!cluster) return std::set<PHG4Hit*>();
  
  if (_do_cache) {
    std::map<RawCluster*,std::set<PHG4Hit*> >::iterator iter =
      _cache_all_truth_hits.find(cluster);
    if (iter != _cache_all_truth_hits.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> truth_hits;
  
  // loop over all the clustered towers
  RawCluster::TowerConstRange begin_end = cluster->get_towers();
  for (RawCluster::TowerConstIterator iter = begin_end.first;
       iter != begin_end.second;
       ++iter) { 

    RawTower* tower = _towers->getTower(iter->first);

    if (_strict) assert(tower);
    else if (!tower) continue;
    
    std::set<PHG4Hit*> new_hits = get_rawtower_eval()->all_truth_hits(tower);

    for (std::set<PHG4Hit*>::iterator iter = new_hits.begin();
	 iter != new_hits.end();
	 ++iter) {

      PHG4Hit* g4hit = *iter;

      if (_strict) assert(g4hit);
      else if (!g4hit) continue;
      
      truth_hits.insert(g4hit);
    }
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(cluster,truth_hits));
  
  return truth_hits;
}
  
std::set<PHG4Particle*> CaloRawClusterEval::all_truth_primaries(RawCluster* cluster) {

  if (_strict) assert(cluster);
  else if (!cluster) return std::set<PHG4Particle*>();
  
  if (_do_cache) {
    std::map<RawCluster*,std::set<PHG4Particle*> >::iterator iter =
      _cache_all_truth_primaries.find(cluster);
    if (iter != _cache_all_truth_primaries.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Particle*> truth_primaries;
  
  // loop over all the clustered towers
  RawCluster::TowerConstRange begin_end = cluster->get_towers();
  for (RawCluster::TowerConstIterator iter = begin_end.first;
       iter != begin_end.second;
       ++iter) {
    
    RawTower* tower = _towers->getTower(iter->first);

    if (_strict) assert(tower);
    else if (!tower) continue;
        
    std::set<PHG4Particle*> new_primaries = _towereval.all_truth_primaries(tower);

    for (std::set<PHG4Particle*>::iterator iter = new_primaries.begin();
	 iter != new_primaries.end();
	 ++iter) {
      PHG4Particle* particle = *iter;

      if (_strict) assert(particle);
      else if (!particle) continue;
      
      truth_primaries.insert(particle);
    }
  }

  if (_do_cache) _cache_all_truth_primaries.insert(make_pair(cluster,truth_primaries));
  
  return truth_primaries;
}

PHG4Particle* CaloRawClusterEval::max_truth_primary_by_energy(RawCluster* cluster) {

  if (_strict) assert(cluster);
  else if (!cluster) return NULL;
  
  if (_do_cache) {
    std::map<RawCluster*,PHG4Particle*>::iterator iter =
      _cache_max_truth_primary_by_energy.find(cluster);
    if (iter != _cache_max_truth_primary_by_energy.end()) {
      return iter->second;
    }
  }
  
  // loop over all primaries associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_primary = NULL;
  float max_e = FLT_MAX*-1.0;
  std::set<PHG4Particle*> primaries = all_truth_primaries(cluster);
  for (std::set<PHG4Particle*>::iterator iter = primaries.begin();
       iter != primaries.end();
       ++iter) {

    PHG4Particle* primary = *iter;

    if (_strict) assert(primary);
    else if (!primary) continue;
    
    float e = get_energy_contribution(cluster,primary);
    if (isnan(e)) continue;
    if (e > max_e) {
      max_e = e;
      max_primary = primary;      
    }
  }

  if (_do_cache) _cache_max_truth_primary_by_energy.insert(make_pair(cluster,max_primary));
  
  return max_primary;
}

std::set<RawCluster*> CaloRawClusterEval::all_clusters_from(PHG4Particle* primary) { 

  if (_strict) assert(primary);
  else if (!primary) return std::set<RawCluster*>();
  
  if (!get_truth_eval()->is_primary(primary)) return std::set<RawCluster*>();

  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict) assert(primary);
  else if (!primary) return std::set<RawCluster*>();
  
  if (_do_cache) {
    std::map<PHG4Particle*,std::set<RawCluster*> >::iterator iter =
      _cache_all_clusters_from_primary.find(primary);
    if (iter != _cache_all_clusters_from_primary.end()) {
      return iter->second;
    }
  }
  
  std::set<RawCluster*> clusters;
  
  // loop over all the clusters
  for (RawClusterContainer::Iterator iter = _clusters->getClusters().first;
       iter != _clusters->getClusters().second;
       ++iter) {

    RawCluster* cluster = iter->second;

    // loop over all truth particles connected to this cluster
    std::set<PHG4Particle*> primaries = all_truth_primaries(cluster);
    for (std::set<PHG4Particle*>::iterator jter = primaries.begin();
	 jter != primaries.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;

      if (_strict) assert(candidate);
      else if (!candidate) continue;
      
      if (candidate->get_track_id() == primary->get_track_id()) {
	clusters.insert(cluster);
      }    
    }
  }

  if (_do_cache) _cache_all_clusters_from_primary.insert(make_pair(primary,clusters));
  
  return clusters;
}

RawCluster* CaloRawClusterEval::best_cluster_from(PHG4Particle* primary) {

  if (_strict) assert(primary);
  else if (!primary) return NULL;
  
  if (!get_truth_eval()->is_primary(primary)) return NULL;

  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict) assert(primary);
  else if (!primary) return NULL;
  
  if (_do_cache) {
    std::map<PHG4Particle*,RawCluster*>::iterator iter =
      _cache_best_cluster_from_primary.find(primary);
    if (iter != _cache_best_cluster_from_primary.end()) {
      return iter->second;
    }
  }
  
  RawCluster* best_cluster = NULL;
  float best_energy = FLT_MAX*-1.0;  
  std::set<RawCluster*> clusters = all_clusters_from(primary);
  for (std::set<RawCluster*>::iterator iter = clusters.begin();
       iter != clusters.end();
       ++iter) {
    RawCluster* cluster = *iter;

    if (_strict) assert(cluster);
    else if (!cluster) continue;
    
    float energy = get_energy_contribution(cluster,primary);
    if(isnan(energy)) continue;
    if (energy > best_energy) {
      best_cluster = cluster;
      best_energy = energy;
    }
  }
 
  if (_do_cache) _cache_best_cluster_from_primary.insert(make_pair(primary,best_cluster));
  
  return best_cluster;
}

// overlap calculations
float CaloRawClusterEval::get_energy_contribution(RawCluster* cluster, PHG4Particle* primary) {

  if (_strict) {
    assert(cluster);
    assert(primary);
  } else if (!cluster||!primary) {
    return NAN;
  }
  
  if (!get_truth_eval()->is_primary(primary)) return NAN;

  // reduce cache misses by using only pointer from PrimaryMap
  primary = get_truth_eval()->get_primary_particle(primary);

  if (_strict) assert(primary);
  else if (!primary) return NULL;
  
  if (_do_cache) {
    std::map<std::pair<RawCluster*,PHG4Particle*>,float>::iterator iter =
      _cache_get_energy_contribution_primary.find(make_pair(cluster,primary));
    if (iter != _cache_get_energy_contribution_primary.end()) {
      return iter->second;
    }
  }

  float energy = 0.0;
  
  std::set<PHG4Particle*> g4particles = all_truth_primaries(cluster);
  if (g4particles.find(primary) != g4particles.end()) {

    std::set<PHG4Hit*> g4hits = all_truth_hits(cluster);
    for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
	 iter != g4hits.end();
	 ++iter) {
      PHG4Hit* g4hit = *iter;
      PHG4Particle* candidate = get_truth_eval()->get_primary_particle(g4hit);

      if (_strict) assert(candidate);
      else if (!candidate) continue;
      
      if (candidate->get_track_id() == primary->get_track_id()) {
	energy += g4hit->get_edep();
      }
    }    
  }
  
  if (_do_cache) _cache_get_energy_contribution_primary.insert(make_pair(make_pair(cluster,primary),energy));
  
  return energy;
}

void CaloRawClusterEval::get_node_pointers(PHCompositeNode* topNode) {

  // need things off of the DST...
  std::string nodename = "CLUSTER_" + _caloname;
  _clusters = findNode::getClass<RawClusterContainer>(topNode,nodename.c_str());
  if (!_clusters) {
    cerr << PHWHERE << " ERROR: Can't find " << nodename << endl;
    exit(-1);
  }

  std::string towername = "TOWER_CALIB_" + _caloname;
  _towers = findNode::getClass<RawTowerContainer>(topNode,towername.c_str());
  if (!_towers) {
    cerr << PHWHERE << " ERROR: Can't find " << towername << endl;
    exit(-1);
  }
  
  return;
}
