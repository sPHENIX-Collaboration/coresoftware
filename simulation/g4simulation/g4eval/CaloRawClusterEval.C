
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

using namespace std;

CaloRawClusterEval::CaloRawClusterEval(PHCompositeNode* topNode, std::string caloname)
  : _topNode(topNode),
    _caloname(caloname),
    _towereval(topNode,caloname),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_primaries(),
    _cache_max_truth_primary_by_energy(),
    _cache_all_clusters_from_primary(),
    _cache_best_cluster_from_primary(),
    _cache_get_energy_contribution_primary() {
}

void CaloRawClusterEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_primaries.clear();
  _cache_max_truth_primary_by_energy.clear();
  _cache_all_clusters_from_primary.clear();
  _cache_best_cluster_from_primary.clear();
  _cache_get_energy_contribution_primary.clear();

  _towereval.next_event(topNode);
  
  _topNode = topNode;  
}

std::set<PHG4Hit*> CaloRawClusterEval::all_truth_hits(RawCluster* cluster) {

  if ((_do_cache) &&
      (_cache_all_truth_hits.find(cluster) != _cache_all_truth_hits.end())) {
    return _cache_all_truth_hits[cluster];
  }
  
  // need things off of the DST...
  std::string towername = "TOWER_" + _caloname;
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(_topNode,towername.c_str());
  if (!towers) {
    cerr << PHWHERE << " ERROR: Can't find " << towers << endl;
    exit(-1);
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the clustered towers
  for (unsigned int itower = 0; itower < cluster->getNTowers(); ++itower) {

    int ieta = cluster->getTowerBin(itower).first;
    int iphi = cluster->getTowerBin(itower).second;
    
    RawTower* tower = towers->getTower(ieta,iphi);
    
    std::set<PHG4Hit*> new_hits = _towereval.all_truth_hits(tower);
    std::set<PHG4Hit*> union_hits;

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   new_hits.begin(),new_hits.end(),
		   std::inserter(union_hits,union_hits.begin()));

    std::swap(truth_hits,union_hits); // swap union into truth_particles    
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(cluster,truth_hits));
  
  return truth_hits;
}
  
std::set<PHG4Particle*> CaloRawClusterEval::all_truth_primaries(RawCluster* cluster) {

  if ((_do_cache) &&
      (_cache_all_truth_primaries.find(cluster) != _cache_all_truth_primaries.end())) {
    return _cache_all_truth_primaries[cluster];
  }

  // need things off of the DST...
  std::string towername = "TOWER_" + _caloname;
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(_topNode,towername.c_str());
  if (!towers) {
    cerr << PHWHERE << " ERROR: Can't find " << towers << endl;
    exit(-1);
  }
  
  std::set<PHG4Particle*> truth_primaries;
  
  // loop over all the clustered towers
  for (unsigned int itower = 0; itower < cluster->getNTowers(); ++itower) {

    int ieta = cluster->getTowerBin(itower).first;
    int iphi = cluster->getTowerBin(itower).second;
    
    RawTower* tower = towers->getTower(ieta,iphi);
    
    std::set<PHG4Particle*> new_primaries = _towereval.all_truth_primaries(tower);
    std::set<PHG4Particle*> union_primaries;

    std::set_union(truth_primaries.begin(),truth_primaries.end(),
		   new_primaries.begin(),new_primaries.end(),
		   std::inserter(union_primaries,union_primaries.begin()));

    std::swap(truth_primaries,union_primaries); // swap union into truth_particles    
  }

  if (_do_cache) _cache_all_truth_primaries.insert(make_pair(cluster,truth_primaries));
  
  return truth_primaries;
}

PHG4Particle* CaloRawClusterEval::max_truth_primary_by_energy(RawCluster* cluster) {

  if ((_do_cache) &&
      (_cache_max_truth_primary_by_energy.find(cluster) != _cache_max_truth_primary_by_energy.end())) {
    return _cache_max_truth_primary_by_energy[cluster];
  }
  
  // loop over all primaries associated with this cluster and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_primary = NULL;
  float max_e = FLT_MIN;
  std::set<PHG4Particle*> primaries = all_truth_primaries(cluster);
  for (std::set<PHG4Particle*>::iterator iter = primaries.begin();
       iter != primaries.end();
       ++iter) {

    PHG4Particle* primary = *iter;
    float e = get_energy_contribution(cluster,primary);
    if (e > max_e) {
      max_e = e;
      max_primary = primary;      
    }
  }

  if (_do_cache) _cache_max_truth_primary_by_energy.insert(make_pair(cluster,max_primary));
  
  return max_primary;
}

std::set<RawCluster*> CaloRawClusterEval::all_clusters_from(PHG4Particle* primary) { 

  CaloTruthEval* trutheval = _towereval.get_truth_eval();
  if (!trutheval->is_primary(primary)) return std::set<RawCluster*>();
  
  if ((_do_cache) &&
      (_cache_all_clusters_from_primary.find(primary) != _cache_all_clusters_from_primary.end())) {
    return _cache_all_clusters_from_primary[primary];
  }
  
  // need things off of the DST...
  RawClusterContainer* clustercontainer = findNode::getClass<RawClusterContainer>(_topNode,"RawClusterContainer");
  if (!clustercontainer) {
    cerr << PHWHERE << " ERROR: Can't find RawClusterContainer" << endl;
    exit(-1);
  }

  std::set<RawCluster*> clusters;
  
  // loop over all the clusters
  for (RawClusterContainer::Iterator iter = clustercontainer->getClusters().first;
       iter != clustercontainer->getClusters().second;
       ++iter) {

    RawCluster* cluster = iter->second;

    // loop over all truth particles connected to this cluster
    std::set<PHG4Particle*> primaries = all_truth_primaries(cluster);
    for (std::set<PHG4Particle*>::iterator jter = primaries.begin();
	 jter != primaries.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;
      if (candidate->get_track_id() == primary->get_track_id()) {
	clusters.insert(cluster);
      }    
    }
  }

  if (_do_cache) _cache_all_clusters_from_primary.insert(make_pair(primary,clusters));
  
  return clusters;
}

RawCluster* CaloRawClusterEval::best_cluster_from(PHG4Particle* primary) {

  CaloTruthEval* trutheval = _towereval.get_truth_eval();
  if (!trutheval->is_primary(primary)) return NULL;
      
  if ((_do_cache) &&
      (_cache_best_cluster_from_primary.find(primary) != _cache_best_cluster_from_primary.end())) {
    return _cache_best_cluster_from_primary[primary];
  }
  
  RawCluster* best_cluster = NULL;
  float best_energy = 0.0;  
  std::set<RawCluster*> clusters = all_clusters_from(primary);
  for (std::set<RawCluster*>::iterator iter = clusters.begin();
       iter != clusters.end();
       ++iter) {
    RawCluster* cluster = *iter;
    float energy = get_energy_contribution(cluster,primary);
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

  CaloTruthEval* trutheval = _towereval.get_truth_eval();
  if (!trutheval->is_primary(primary)) return NAN;
  
  if ((_do_cache) &&
      (_cache_get_energy_contribution_primary.find(make_pair(cluster,primary)) !=
       _cache_get_energy_contribution_primary.end())) {
    return _cache_get_energy_contribution_primary[make_pair(cluster,primary)];
  }
  
  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(cluster);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* g4hit = *iter;
    if (g4hit->get_trkid() == primary->get_track_id()) {
      energy += g4hit->get_edep();
    }
  }

  if (_do_cache) _cache_get_energy_contribution_primary.insert(make_pair(make_pair(cluster,primary),energy));
  
  return energy;
}
