
#include "CaloTruthEval.h"

#include "BaseTruthEval.h"

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <cstdlib>
#include <set>
#include <map>
#include <float.h>
#include <iostream>
#include <cassert>

using namespace std;

CaloTruthEval::CaloTruthEval(PHCompositeNode* topNode,std::string caloname)
  : _basetrutheval(topNode),
    _caloname(caloname),
    _truthinfo(NULL),
    _g4hits(NULL),
    _strict(false),
    _verbosity(1),
    _errors(0),
    _do_cache(true),
    _cache_all_truth_hits_g4particle(),
    _cache_get_primary_particle_g4hit(),
    _cache_get_shower_from_primary(),
    _cache_get_shower_moliere_radius(),
    _cache_get_shower_energy_deposit() {
  get_node_pointers(topNode);
}

CaloTruthEval::~CaloTruthEval() {
  if (_verbosity > 0) {
    if ((_errors > 0)||(_verbosity > 1)) {
      cout << "CaloTruthEval::~CaloTruthEval() - Error Count: " << _errors << endl;
    }
  }
}

void CaloTruthEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits_g4particle.clear();
  _cache_get_primary_particle_g4hit.clear();
  _cache_get_shower_from_primary.clear();
  _cache_get_shower_moliere_radius.clear();
  _cache_get_shower_energy_deposit.clear();

  _basetrutheval.next_event(topNode);
  
  get_node_pointers(topNode);
}

std::set<PHG4Hit*> CaloTruthEval::all_truth_hits(PHG4Particle* particle) {

  if (!has_node_pointers()) {++_errors; return std::set<PHG4Hit*>();}
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return std::set<PHG4Hit*>();}
  
  if (_do_cache) {
    std::map<PHG4Particle*,std::set<PHG4Hit*> >::iterator iter =
      _cache_all_truth_hits_g4particle.find(particle);
    if (iter != _cache_all_truth_hits_g4particle.end()) {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits
  for (PHG4HitContainer::ConstIterator g4iter = _g4hits->getHits().first;
       g4iter != _g4hits->getHits().second;
       ++g4iter) {

    PHG4Hit* g4hit = g4iter->second;
    if (is_g4hit_from_particle(g4hit,particle)) continue;
    truth_hits.insert(g4hit);
  }
  
  if (_do_cache) _cache_all_truth_hits_g4particle.insert(make_pair(particle,truth_hits));
  
  return truth_hits;
}

PHG4Particle* CaloTruthEval::get_parent_particle(PHG4Hit* g4hit) {
  return _basetrutheval.get_particle(g4hit);
}

PHG4Particle* CaloTruthEval::get_primary_particle(PHG4Particle* particle) {
  return _basetrutheval.get_primary(particle);
}

PHG4Particle* CaloTruthEval::get_primary_particle(PHG4Hit* g4hit) {

  if (!has_node_pointers()) {++_errors; return NULL;}
  
  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; return NULL;}
  
  if (_do_cache) {
    std::map<PHG4Hit*,PHG4Particle*>::iterator iter =
      _cache_get_primary_particle_g4hit.find(g4hit);
    if (iter != _cache_get_primary_particle_g4hit.end()) {
      return iter->second;
    }
  }
  
  PHG4Particle* primary = _basetrutheval.get_primary(g4hit);

  if (_do_cache) _cache_get_primary_particle_g4hit.insert(make_pair(g4hit,primary));

  if (_strict) {assert(primary);}
  else if (!primary) {++_errors;}
  
  return primary;
}

int CaloTruthEval::get_embed(PHG4Particle* particle) {
  return _basetrutheval.get_embed(particle);
}

PHG4VtxPoint* CaloTruthEval::get_vertex(PHG4Particle* particle) {
  return _basetrutheval.get_vertex(particle);
}

bool CaloTruthEval::is_primary(PHG4Particle* particle) {
  return _basetrutheval.is_primary(particle);
}

std::set<PHG4Hit*> CaloTruthEval::get_shower_from_primary(PHG4Particle* primary) {

  if (!has_node_pointers()) {++_errors; return std::set<PHG4Hit*>();}
  
  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return std::set<PHG4Hit*>();}
  
  if (!is_primary(primary)) return std::set<PHG4Hit*>();

  primary = get_primary_particle(primary);

  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return std::set<PHG4Hit*>();}
  
  if (_do_cache) {
    std::map<PHG4Particle*,std::set<PHG4Hit*> >::iterator iter =
      _cache_get_shower_from_primary.find(primary);
    if (iter != _cache_get_shower_from_primary.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits
  for (PHG4HitContainer::ConstIterator g4iter = _g4hits->getHits().first;
       g4iter != _g4hits->getHits().second;
       ++g4iter) {

    PHG4Hit* g4hit = g4iter->second;
    PHG4Particle* candidate = get_primary_particle(g4hit);

    if (_strict) {assert(candidate);}
    else if (!candidate) {++_errors; continue;}

    if (!are_same_particle(candidate,primary)) continue;
    truth_hits.insert(g4hit);
  }

  if (_do_cache) _cache_get_shower_from_primary.insert(make_pair(primary,truth_hits));
  
  return truth_hits;
}

// moliere (90% containment) radius of scintilator hits
float CaloTruthEval::get_shower_moliere_radius(PHG4Particle* primary) {

  if (!has_node_pointers()) {++_errors; return NAN;}
  
  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return NAN;}
  
  if (!is_primary(primary)) return NAN;

  primary = get_primary_particle(primary);

  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return NAN;}
  
  if (_do_cache) {
    std::map<PHG4Particle*,float>::iterator iter =
      _cache_get_shower_moliere_radius.find(primary);
    if (iter != _cache_get_shower_moliere_radius.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> g4hits = get_shower_from_primary(primary);

  std::multimap<float,float> radii_energy_mmap;
  float shower_e = 0.0;
  
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {

    PHG4Hit *g4hit = (*iter);
	      
    // momentum vector
    /// \todo for charged particles remove magnetic field bend
    float p_x = primary->get_px();
    float p_y = primary->get_py();
    float p_z = primary->get_pz();
    float p   = sqrt(pow(p_x,2)+pow(p_y,2)+pow(p_z,2));

    // relative position vector (vertex-to-ghit)
    PHG4VtxPoint* vtx = get_vertex(primary);
    float d_x = vtx->get_x() - g4hit->get_avg_x();
    float d_y = vtx->get_y() - g4hit->get_avg_y();
    float d_z = vtx->get_z() - g4hit->get_avg_z();
    float d   = sqrt(pow(d_x,2)+pow(d_y,2)+pow(d_z,2));
    
    // angle between them
    float phi = acos( (p_x*d_x+p_y*d_y+p_z*d_z)/p/d );

    // distance between them at ghit
    float r = d*sin(phi); 
    float edep = g4hit->get_edep();
    shower_e += edep;
    
    radii_energy_mmap.insert(make_pair(r,edep));
  }

  float sum_e = 0.0;
  float frac_e = 0.0;

  float r_in = 0.0;
  float r_out = 0.0;

  for(std::multimap<float,float>::iterator iter = radii_energy_mmap.begin();
      iter != radii_energy_mmap.end();
      iter++) {
    r_out = iter->first;
    sum_e = sum_e + iter->second;
    frac_e = sum_e / shower_e;

    if (frac_e > 0.90) break;
    
    r_in = r_out;
  }

  float radius = 0.5*(r_in+r_out);

  if (_do_cache) _cache_get_shower_moliere_radius.insert(make_pair(primary,radius));
  
  return radius;
}

float CaloTruthEval::get_shower_energy_deposit(PHG4Particle* primary) {

  if (!has_node_pointers()) {++_errors; return NAN;}
  
  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return NAN;}
  
  if (!is_primary(primary)) return NAN;

  primary = get_primary_particle(primary);

  if (_strict) {assert(primary);}
  else if (!primary) {++_errors; return NAN;}
  
  if (_do_cache) {
    std::map<PHG4Particle*,float>::iterator iter =
      _cache_get_shower_energy_deposit.find(primary);
    if (iter != _cache_get_shower_energy_deposit.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> g4hits = get_shower_from_primary(primary);

  float shower_e = 0.0;  
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* g4hit = (*iter);
    shower_e += g4hit->get_edep();
  }

  if (_do_cache) _cache_get_shower_energy_deposit.insert(make_pair(primary,shower_e));
  
  return shower_e;
}

bool CaloTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle) {
  return _basetrutheval.is_g4hit_from_particle(g4hit,particle);
}

bool CaloTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2) {
  return _basetrutheval.are_same_particle(p1,p2);
}

bool CaloTruthEval::are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2) {
  return _basetrutheval.are_same_vertex(vtx1,vtx2);
}

void CaloTruthEval::get_node_pointers(PHCompositeNode *topNode) {

  // need things off of the DST...
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
 
  std::string name = "G4HIT_" + _caloname;
  _g4hits = findNode::getClass<PHG4HitContainer>(topNode,name.c_str());
  
  return;
}

bool CaloTruthEval::has_node_pointers() {

  if (_strict) assert(_truthinfo);
  else if (!_truthinfo) return false;

  if (_strict) assert(_g4hits);
  if (!_g4hits) return false;

  return true;
}
