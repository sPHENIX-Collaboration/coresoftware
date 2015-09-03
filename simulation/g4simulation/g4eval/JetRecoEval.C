
#include "JetRecoEval.h"

#include "JetTruthEval.h"

#include "SvtxEvalStack.h"
#include "CaloEvalStack.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>

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

JetRecoEval::JetRecoEval(PHCompositeNode* topNode,
			 std::string recojetname,
			 std::string truthjetname)
  : _jettrutheval(topNode,recojetname,truthjetname),
    _recojetname(recojetname),
    _truthjetname(truthjetname),
    _recojets(NULL),
    _truthjets(NULL),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_jets(),
    _cache_max_truth_jet_by_energy(),
    _cache_all_truth_particles(),
    _cache_all_jets_from(),
    _cache_best_jet_from(),
    _cache_get_energy_contribution() {
  get_node_pointers(topNode);
}

void JetRecoEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_jets.clear();
  _cache_max_truth_jet_by_energy.clear();
  _cache_all_truth_particles.clear();
  _cache_all_jets_from.clear();
  _cache_best_jet_from.clear();
  _cache_get_energy_contribution.clear();

  _jettrutheval.next_event(topNode);
  
  get_node_pointers(topNode);
}

std::set<PHG4Hit*> JetRecoEval::all_truth_hits(Jet* recojet) {

  if (_do_cache) {
    std::map<Jet*,std::set<PHG4Hit*> >::iterator iter =
      _cache_all_truth_hits.find(recojet);
    if (iter != _cache_all_truth_hits.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> truth_hits;

  // loop over all the clustered towers
  for (unsigned int itower = 0; itower < cluster->getNTowers(); ++itower) {

    int ieta = cluster->getTowerBin(itower).first;
    int iphi = cluster->getTowerBin(itower).second;
    
    RawTower* tower = _towers->getTower(ieta,iphi);
    
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
  
std::set<PHG4Particle*> JetRecoEval::all_truth_primaries(RawCluster* cluster) {

  if (_do_cache) {
    std::map<RawCluster*,std::set<PHG4Particle*> >::iterator iter =
      _cache_all_truth_primaries.find(cluster);
    if (iter != _cache_all_truth_primaries.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Particle*> truth_primaries;
  
  // loop over all the clustered towers
  for (unsigned int itower = 0; itower < cluster->getNTowers(); ++itower) {

    int ieta = cluster->getTowerBin(itower).first;
    int iphi = cluster->getTowerBin(itower).second;
    
    RawTower* tower = _towers->getTower(ieta,iphi);
    
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

PHG4Particle* JetRecoEval::max_truth_primary_by_energy(RawCluster* cluster) {

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
    float e = get_energy_contribution(cluster,primary);
    if (e > max_e) {
      max_e = e;
      max_primary = primary;      
    }
  }

  if (_do_cache) _cache_max_truth_primary_by_energy.insert(make_pair(cluster,max_primary));
  
  return max_primary;
}

std::set<RawCluster*> JetRecoEval::all_clusters_from(PHG4Particle* primary) { 

  CaloTruthEval* trutheval = _towereval.get_truth_eval();
  if (!trutheval->is_primary(primary)) return std::set<RawCluster*>();
  
  if ((_do_cache) &&
      (_cache_all_clusters_from_primary.find(primary) != _cache_all_clusters_from_primary.end())) {
    return _cache_all_clusters_from_primary[primary];
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
      if (candidate->get_track_id() == primary->get_track_id()) {
	clusters.insert(cluster);
      }    
    }
  }

  if (_do_cache) _cache_all_clusters_from_primary.insert(make_pair(primary,clusters));
  
  return clusters;
}

RawCluster* JetRecoEval::best_cluster_from(PHG4Particle* primary) {

  CaloTruthEval* trutheval = _towereval.get_truth_eval();
  if (!trutheval->is_primary(primary)) return NULL;
      
  if (_do_cache) {
    std::map<PHG4Particle*,RawCluster*>::iterator iter =
      _cache_best_cluster_from_primary.find(primary);
    if (iter != _cache_best_cluster_from_primary.end()) {
      return iter->second;
    }
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
float JetRecoEval::get_energy_contribution(RawCluster* cluster, PHG4Particle* primary) {

  CaloTruthEval* trutheval = _towereval.get_truth_eval();
  if (!trutheval->is_primary(primary)) return NAN;
  
  if (_do_cache) {
    std::map<std::pair<RawCluster*,PHG4Particle*>,float>::iterator iter =
      _cache_get_energy_contribution_primary.find(make_pair(cluster,primary));
    if (iter != _cache_get_energy_contribution_primary.end()) {
      return iter->second;
    }
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

void JetRecoEval::get_node_pointers(PHCompositeNode* topNode) {

  // need things off of the DST...
  _recojets = findNode::getClass<JetMap>(topNode,_reconamejet.c_str());
  if (!_recojets) {
    cerr << PHWHERE << " ERROR: Can't find " << _recojetname << endl;
    exit(-1);
  }

  _truthjets = findNode::getClass<JetMap>(topNode,_truthnamejet.c_str());
  if (!_truthjets) {
    cerr << PHWHERE << " ERROR: Can't find " << _truthjetname << endl;
    exit(-1);
  }

  return;
}
