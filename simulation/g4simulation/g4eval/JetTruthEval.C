
#include "JetTruthEval.h"

#include <phool/getClass.h>
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
#include <cassert>

using namespace std;

JetTruthEval::JetTruthEval(PHCompositeNode* topNode,
			   std::string truthjetname)
  : _truthjetname(truthjetname),
    _svtxevalstack(topNode),   
    _cemcevalstack(topNode,"CEMC"),
    _hcalinevalstack(topNode,"HCALIN"),
    _hcaloutevalstack(topNode,"HCALOUT"),
    _femcevalstack(topNode,"FEMC"),
    _fhcalevalstack(topNode,"FHCAL"),
    _truthinfo(NULL),
    _truthjets(NULL),
    _strict(false),
    _verbosity(1),
    _errors(0),
    _do_cache(true),
    _cache_all_truth_particles(),
    _cache_all_truth_hits(),
    _cache_get_truth_jet() {
  get_node_pointers(topNode);
}

JetTruthEval::~JetTruthEval() {
  if (_verbosity > 0) {
    if ((_errors > 0)||(_verbosity > 1)) {
      cout << "JetTruthEval::~JetTruthEval() - Error Count: " << _errors << endl;
    }
  }
}

void JetTruthEval::next_event(PHCompositeNode* topNode) {

  _svtxevalstack.next_event(topNode);
  _cemcevalstack.next_event(topNode);
  _hcalinevalstack.next_event(topNode);
  _hcaloutevalstack.next_event(topNode);
  _femcevalstack.next_event(topNode);
  _fhcalevalstack.next_event(topNode);
  
  _cache_all_truth_particles.clear();
  _cache_all_truth_hits.clear();
  _cache_get_truth_jet.clear();

  get_node_pointers(topNode);
}

std::set<PHG4Particle*> JetTruthEval::all_truth_particles(Jet* truthjet) {

  if (_strict) {assert(truthjet);}
  else if (!truthjet) {++_errors; return std::set<PHG4Particle*>();}
  
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

    PHG4Particle* truth_particle = _truthinfo->GetParticle(index);
    
    if (_strict) {assert(truth_particle);}
    else if (!truth_particle) {++_errors; continue;}
    
    truth_particles.insert(truth_particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(truthjet,truth_particles));

  return truth_particles;
}

std::set<PHG4Hit*> JetTruthEval::all_truth_hits(Jet* truthjet) {

  if (_strict) {assert(truthjet);}
  else if (!truthjet) {++_errors; return std::set<PHG4Hit*>();}
  
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

    if (_strict) {assert(particle);}
    else if (!particle) {++_errors; continue;}
    
    // ask the svtx truth eval to backtrack the particles to g4hits
    SvtxTruthEval* svtx_truth_eval = _svtxevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> svtx_g4hits = svtx_truth_eval->all_truth_hits(particle);

    for (std::set<PHG4Hit*>::iterator jter = svtx_g4hits.begin();
	 jter != svtx_g4hits.end();
	 ++jter) {

      PHG4Hit* g4hit = *jter;

      if (_strict) {assert(g4hit);}
      else if (!g4hit) {++_errors; continue;}
      
      truth_hits.insert(g4hit);
    }
    
    // ask the cemc truth eval to backtrack the primary to g4hits
    CaloTruthEval* cemc_truth_eval = _cemcevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> cemc_g4hits = cemc_truth_eval->get_shower_hits_from_primary(particle);

    for (std::set<PHG4Hit*>::iterator jter = cemc_g4hits.begin();
	 jter != cemc_g4hits.end();
	 ++jter) {
      
      PHG4Hit* g4hit = *jter;

      if (_strict) {assert(g4hit);}
      else if (!g4hit) {++_errors; continue;}
      
      truth_hits.insert(g4hit);
    }
    
    // ask the hcalin truth eval to backtrack the primary to g4hits
    CaloTruthEval* hcalin_truth_eval = _hcalinevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> hcalin_g4hits = hcalin_truth_eval->get_shower_hits_from_primary(particle);

    for (std::set<PHG4Hit*>::iterator jter = hcalin_g4hits.begin();
	 jter != hcalin_g4hits.end();
	 ++jter) {

      PHG4Hit* g4hit = *jter;

      if (_strict) {assert(g4hit);}
      else if (!g4hit) {++_errors; continue;}
      
      truth_hits.insert(g4hit);
    }
    
    // ask the hcalout truth eval to backtrack the primary to g4hits
    CaloTruthEval* hcalout_truth_eval = _hcaloutevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> hcalout_g4hits = hcalout_truth_eval->get_shower_hits_from_primary(particle);

    for (std::set<PHG4Hit*>::iterator jter = hcalout_g4hits.begin();
	 jter != hcalout_g4hits.end();
	 ++jter) {

      PHG4Hit* g4hit = *jter;

      if (_strict) {assert(g4hit);}
      else if (!g4hit) {++_errors; continue;}
      
      truth_hits.insert(g4hit);      
    }    

    // ask the femc truth eval to backtrack the primary to g4hits
    CaloTruthEval* femc_truth_eval = _femcevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> femc_g4hits = femc_truth_eval->get_shower_hits_from_primary(particle);

    for (std::set<PHG4Hit*>::iterator jter = femc_g4hits.begin();
	 jter != femc_g4hits.end();
	 ++jter) {

      PHG4Hit* g4hit = *jter;

      if (_strict) {assert(g4hit);}
      else if (!g4hit) {++_errors; continue;}
      
      truth_hits.insert(g4hit);      
    }    

    // ask the fhcal truth eval to backtrack the primary to g4hits
    CaloTruthEval* fhcal_truth_eval = _fhcalevalstack.get_truth_eval(); 
    std::set<PHG4Hit*> fhcal_g4hits = fhcal_truth_eval->get_shower_hits_from_primary(particle);

    for (std::set<PHG4Hit*>::iterator jter = fhcal_g4hits.begin();
	 jter != fhcal_g4hits.end();
	 ++jter) {

      PHG4Hit* g4hit = *jter;

      if (_strict) {assert(g4hit);}
      else if (!g4hit) {++_errors; continue;}
      
      truth_hits.insert(g4hit);      
    }    

  }
  
  if (_do_cache) _cache_all_truth_hits.insert(make_pair(truthjet,truth_hits));
  
  return truth_hits;
}

Jet* JetTruthEval::get_truth_jet(PHG4Particle* particle) {

  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return NULL;}
  
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
      
      PHG4Particle* constituent = _truthinfo->GetParticle( index );
      if (_strict) {assert(constituent);}
      else if (!constituent) {++_errors; continue;}

      if (get_svtx_eval_stack()->get_truth_eval()->are_same_particle(constituent,particle)) {
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
