
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

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

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
  : _jettrutheval(topNode,truthjetname),
    _recojetname(recojetname),
    _truthjetname(truthjetname),
    _recojets(NULL),
    _truthjets(NULL),
    _trackmap(NULL),
    _cemctowers(NULL),
    _cemcclusters(NULL),
    _hcalintowers(NULL),
    _hcalinclusters(NULL),
    _hcalouttowers(NULL),
    _hcaloutclusters(NULL),
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

  // loop over all the jet constituents, backtrack each reco object to the
  // truth hits and combine with other consituents

  for (Jet::ConstIter iter = recojet->begin_comp();
       iter != recojet->end_comp();
       ++iter) {
    Jet::SRC source = iter->first;
    unsigned int index = iter->second;

    std::set<PHG4Hit*> new_hits;
    
    if (source == Jet::TRACK) {
      SvtxTrack* track = _trackmap->get(index);     
      new_hits = get_svtx_eval_stack()->get_track_eval()->all_truth_hits(track);      
    } else if (source == Jet::CEMC_TOWER) {
      RawTower* tower = _cemctowers->getTower(index);
      new_hits = get_cemc_eval_stack()->get_rawtower_eval()->all_truth_hits(tower);      
    } else if (source == Jet::CEMC_CLUSTER) {
      RawCluster* cluster = _cemcclusters->getCluster(index);
      new_hits = get_cemc_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster);      
    } else if (source == Jet::HCALIN_TOWER) {
      RawTower* tower = _hcalintowers->getTower(index);
      new_hits = get_hcalin_eval_stack()->get_rawtower_eval()->all_truth_hits(tower); 
    } else if (source == Jet::HCALIN_CLUSTER) {
      RawCluster* cluster = _hcalinclusters->getCluster(index);
      new_hits = get_hcalin_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster); 
    } else if (source == Jet::HCALOUT_TOWER) {
      RawTower* tower = _hcalouttowers->getTower(index);
      new_hits = get_hcalout_eval_stack()->get_rawtower_eval()->all_truth_hits(tower); 
    } else if (source == Jet::HCALOUT_CLUSTER) {
      RawCluster* cluster = _hcaloutclusters->getCluster(index);
      new_hits = get_hcalout_eval_stack()->get_rawcluster_eval()->all_truth_hits(cluster); 
    }
    
    std::set<PHG4Hit*> union_hits;

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   new_hits.begin(),new_hits.end(),
		   std::inserter(union_hits,union_hits.begin()));

    std::swap(truth_hits,union_hits); // swap union into truth_particles    
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(recojet,truth_hits));
  
  return truth_hits;
}

std::set<PHG4Particle*> JetRecoEval::all_truth_particles(Jet* recojet) {
 
  if (_do_cache) {
    std::map<Jet*,std::set<PHG4Particle*> >::iterator iter =
      _cache_all_truth_particles.find(recojet);
    if (iter != _cache_all_truth_particles.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Particle*> truth_particles;

  // loop over all the jet constituents, backtrack each reco object to the
  // truth hits and combine with other consituents

  for (Jet::ConstIter iter = recojet->begin_comp();
       iter != recojet->end_comp();
       ++iter) {
    Jet::SRC source = iter->first;
    unsigned int index = iter->second;

    std::set<PHG4Particle*> new_particles;
    
    if (source == Jet::TRACK) {
      SvtxTrack* track = _trackmap->get(index);     
      new_particles = get_svtx_eval_stack()->get_track_eval()->all_truth_particles(track);      
    } else if (source == Jet::CEMC_TOWER) {
      RawTower* tower = _cemctowers->getTower(index);
      new_particles = get_cemc_eval_stack()->get_rawtower_eval()->all_truth_primaries(tower);      
    } else if (source == Jet::CEMC_CLUSTER) {
      RawCluster* cluster = _cemcclusters->getCluster(index);
      new_particles = get_cemc_eval_stack()->get_rawcluster_eval()->all_truth_primaries(cluster);      
    } else if (source == Jet::HCALIN_TOWER) {
      RawTower* tower = _hcalintowers->getTower(index);
      new_particles = get_hcalin_eval_stack()->get_rawtower_eval()->all_truth_primaries(tower); 
    } else if (source == Jet::HCALIN_CLUSTER) {
      RawCluster* cluster = _hcalinclusters->getCluster(index);
      new_particles = get_hcalin_eval_stack()->get_rawcluster_eval()->all_truth_primaries(cluster); 
    } else if (source == Jet::HCALOUT_TOWER) {
      RawTower* tower = _hcalouttowers->getTower(index);
      new_particles = get_hcalout_eval_stack()->get_rawtower_eval()->all_truth_primaries(tower); 
    } else if (source == Jet::HCALOUT_CLUSTER) {
      RawCluster* cluster = _hcaloutclusters->getCluster(index);
      new_particles = get_hcalout_eval_stack()->get_rawcluster_eval()->all_truth_primaries(cluster); 
    }
    
    std::set<PHG4Particle*> union_hits;

    std::set_union(truth_particles.begin(),truth_particles.end(),
		   new_particles.begin(),new_particles.end(),
		   std::inserter(union_hits,union_hits.begin()));

    std::swap(truth_particles,union_hits); // swap union into truth_particles    
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(recojet,truth_particles));
  
  return truth_particles;
}

std::set<Jet*> JetRecoEval::all_truth_jets(Jet* recojet) {

  if (_do_cache) {
    std::map<Jet*,std::set<Jet*> >::iterator iter =
      _cache_all_truth_jets.find(recojet);
    if (iter != _cache_all_truth_jets.end()) {
      return iter->second;
    }
  }
  
  std::set<Jet*> truth_jets;
  
  // get all truth particles...
  std::set<PHG4Particle*> particles = all_truth_particles(recojet);
  
  // backtrack from the truth particles to the truth jets...
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter) {
    PHG4Particle* particle = *iter;

    Jet* truth_jet = _jettrutheval.get_truth_jet(particle);
    truth_jets.insert(truth_jet);
  }

  if (_do_cache) _cache_all_truth_jets.insert(make_pair(recojet,truth_jets));
  
  return truth_jets;
}

Jet* JetRecoEval::max_truth_jet_by_energy(Jet* recojet) {
  return NULL;
}

std::set<Jet*> JetRecoEval::all_jets_from(Jet* truthjet) {
  return std::set<Jet*>();
}

Jet* JetRecoEval::best_jet_from(Jet* truthjet) {
  return NULL;
}

// overlap calculations
float JetRecoEval::get_energy_contribution(Jet* recojet, Jet* truthjet) {
  return 0.0;
}

void JetRecoEval::get_node_pointers(PHCompositeNode* topNode) {

  // need things off of the DST...
  _recojets = findNode::getClass<JetMap>(topNode,_recojetname.c_str());
  if (!_recojets) {
    cerr << PHWHERE << " ERROR: Can't find " << _recojetname << endl;
    exit(-1);
  }

  _truthjets = findNode::getClass<JetMap>(topNode,_truthjetname.c_str());
  if (!_truthjets) {
    cerr << PHWHERE << " ERROR: Can't find " << _truthjetname << endl;
    exit(-1);
  }

  _trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  _cemctowers = findNode::getClass<RawTowerContainer>(topNode,"TOWERS_CEMC");
  _hcalintowers = findNode::getClass<RawTowerContainer>(topNode,"TOWERS_HCALIN");
  _hcalouttowers = findNode::getClass<RawTowerContainer>(topNode,"TOWERS_HCALOUT");
  _cemcclusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTERS_CEMC");
  _hcalinclusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTERS_HCALIN");
  _hcaloutclusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTERS_HCALOUT");

  return;
}
