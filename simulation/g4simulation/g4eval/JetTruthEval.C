
#include "JetTruthEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <cstdlib>
#include <map>
#include <float.h>
#include <algorithm>

using namespace std;

JetTruthEval::JetTruthEval(PHCompositeNode* topNode,
			   std::string truthjetname)
  : _truthjetname(truthjetname),
    _svtxevalstack(topNode),   
    _cemcevalstack(topNode,"CEMC"),
    _hcalinevalstack(topNode,"HCALIN"),
    _hcaloutevalstack(topNode,"HCALOUT"),
    _truthinfo(NULL),
    _truthjets(NULL),
    _do_cache(true),
    _cache_all_truth_particles(),
    _cache_all_truth_hits(),
    _cache_get_truth_jet() {
  get_node_pointers(topNode);
}

void JetTruthEval::next_event(PHCompositeNode* topNode) {

  _svtxevalstack.next_event(topNode);
  _cemcevalstack.next_event(topNode);
  _hcalinevalstack.next_event(topNode);
  _hcaloutevalstack.next_event(topNode);
  
  _cache_all_truth_particles.clear();
  _cache_all_truth_hits.clear();
  _cache_get_truth_jet.clear();

  get_node_pointers(topNode);
}

std::set<PHG4Particle*> JetTruthEval::all_truth_particles(Jet* truthjet) {

  if (_do_cache) {
    std::map<Jet*,std::set<PHG4Particle*> >::iterator iter =
      _cache_all_truth_particles.find(truthjet);
    if (iter != _cache_all_truth_particles.end()) {
      return iter->second;
    }
  }
    
  std::set<PHG4Particle*> truth_particles;

  // loop over all the entries in the truthjet
  for (Jet::ConstIter iter = truthjet->begin_comp();
       iter != truthjet->end_comp();
       ++iter) {
    Jet::SRC source = iter->first;
    unsigned int index = iter->second;
    if (source != Jet::PARTICLE) {
      cout << PHWHERE << " truth jet contains something other than particles!" << endl;
      exit(-1);
    }

    PHG4Particle* truth_particle = _truthinfo->GetHit(index);
    truth_particles.insert(truth_particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(truthjet,truth_particles));

  return truth_particles;
}

std::set<PHG4Hit*> JetTruthEval::all_truth_hits(Jet* truthjet) {

  if (_do_cache) {
    std::map<Jet*,std::set<PHG4Hit*> >::iterator iter =
      _cache_all_truth_hits.find(truthjet);
    if (iter != _cache_all_truth_hits.end()) {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  std::set<PHG4Particle*> truth_particles = all_truth_particles(truthjet);
  
  // loop over all the entries in the truthjet
  for (std::set<PHG4Particle*>::iterator iter = truth_particles.begin();
       iter != truth_particles.end();
       ++iter) {
    PHG4Particle* particle = *iter;

    // ask the svtx truth eval to backtrack the particles to g4hits
    SvtxTruthEval* svtx_truth_eval = _svtxevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> svtx_g4hits = svtx_truth_eval->all_truth_hits(particle);

    std::set<PHG4Hit*> union_g4hits;
    
    std::set_union(truth_hits.begin(),truth_hits.end(),
		   svtx_g4hits.begin(),svtx_g4hits.end(),
		   std::inserter(union_g4hits,union_g4hits.begin()));
    
    std::swap(truth_hits,union_g4hits); // swap union into truth_hits
    union_g4hits.clear();
    
    // ask the cemc truth eval to backtrack the primary to g4hits
    CaloTruthEval* cemc_truth_eval = _cemcevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> cemc_g4hits = cemc_truth_eval->get_shower_from_primary(particle);

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   cemc_g4hits.begin(),cemc_g4hits.end(),
		   std::inserter(union_g4hits,union_g4hits.begin()));
    
    std::swap(truth_hits,union_g4hits); // swap union into truth_hits
    union_g4hits.clear();
    
    // ask the hcalin truth eval to backtrack the primary to g4hits
    CaloTruthEval* hcalin_truth_eval = _hcalinevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> hcalin_g4hits = hcalin_truth_eval->get_shower_from_primary(particle);

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   hcalin_g4hits.begin(),hcalin_g4hits.end(),
		   std::inserter(union_g4hits,union_g4hits.begin()));
    
    std::swap(truth_hits,union_g4hits); // swap union into truth_hits
    union_g4hits.clear();
    
    // ask the hcalout truth eval to backtrack the primary to g4hits
    CaloTruthEval* hcalout_truth_eval = _hcaloutevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> hcalout_g4hits = hcalout_truth_eval->get_shower_from_primary(particle);

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   hcalout_g4hits.begin(),hcalout_g4hits.end(),
		   std::inserter(union_g4hits,union_g4hits.begin()));
    
    std::swap(truth_hits,union_g4hits); // swap union into truth_hits
    union_g4hits.clear();    
  }
  
  if (_do_cache) _cache_all_truth_hits.insert(make_pair(truthjet,truth_hits));
  
  return truth_hits;
}

Jet* JetTruthEval::get_truth_jet(PHG4Particle* particle) {

  if (_do_cache) {
    std::map<PHG4Particle*,Jet*>::iterator iter =
      _cache_get_truth_jet.find(particle);
    if (iter != _cache_get_truth_jet.end()) {
      return iter->second;
    }
  }

  Jet* truth_jet = NULL;
  
  // loop over all jets and look for this particle...
  for (JetMap::Iter iter = _truthjets->begin();
       iter != _truthjets->end();
       ++iter) {
    Jet* candidate = iter->second;

    // loop over all consituents and look for this particle    
    for (Jet::ConstIter jter = candidate->begin_comp();
	 jter != candidate->end_comp();
	 ++jter) {
      unsigned int index = jter->second;
      if (index == (unsigned int)particle->get_track_id()) {
	truth_jet = candidate;
	break;
      }
    }
    
    if (truth_jet) break;
  }  
  
  if (_do_cache) _cache_get_truth_jet.insert(make_pair(particle,truth_jet));
  
  return truth_jet;
}

void JetTruthEval::get_node_pointers(PHCompositeNode* topNode) {

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  _truthjets = findNode::getClass<JetMap>(topNode,_truthjetname.c_str());
  if (!_truthjets) {
    cerr << PHWHERE << " ERROR: Can't find " << _truthjetname << endl;
    exit(-1);
  }

  return;
}
