
#include "CaloTruthEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <cstdlib>
#include <set>
#include <map>
#include <float.h>

using namespace std;

CaloTruthEval::CaloTruthEval(PHCompositeNode* topNode,std::string caloname)
  : _caloname(caloname),
    _truthinfo(NULL),
    _g4hits(NULL),
    _do_cache(true),
    _cache_all_truth_hits_g4particle(),
    _cache_is_primary(),
    _cache_get_shower_from_primary(),
    _cache_get_shower_moliere_radius(),
    _cache_get_shower_energy_deposit() {
  get_node_pointers(topNode);
}

void CaloTruthEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits_g4particle.clear();
  _cache_is_primary.clear();
  _cache_get_shower_from_primary.clear();
  _cache_get_shower_moliere_radius.clear();
  _cache_get_shower_energy_deposit.clear();
  
  get_node_pointers(topNode);
}

std::set<PHG4Hit*> CaloTruthEval::all_truth_hits(PHG4Particle* particle) {

  if ((_do_cache) &&
      (_cache_all_truth_hits_g4particle.find(particle) != _cache_all_truth_hits_g4particle.end())) {
    return _cache_all_truth_hits_g4particle[particle];
  }

  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits
  for (PHG4HitContainer::ConstIterator g4iter = _g4hits->getHits().first;
       g4iter != _g4hits->getHits().second;
       ++g4iter) {

    PHG4Hit* g4hit = g4iter->second;
    if (g4hit->get_trkid() != particle->get_track_id()) continue;
    truth_hits.insert(g4hit);
  }
  
  if (_do_cache) _cache_all_truth_hits_g4particle.insert(make_pair(particle,truth_hits));
  
  return truth_hits;
}

PHG4Particle* CaloTruthEval::get_parent_particle(PHG4Hit* g4hit) {

  PHG4Particle* particle = _truthinfo->GetHit( g4hit->get_trkid() );
  return particle;
}

PHG4Particle* CaloTruthEval::get_primary_particle(PHG4Particle* particle) {

  PHG4Particle* returnval = particle;
  if (returnval->get_primary_id() != (int)(0xFFFFFFFF)) {
    returnval = _truthinfo->GetHit( particle->get_primary_id() );
  }
  
  return returnval;
}


PHG4Particle* CaloTruthEval::get_primary_particle(PHG4Hit* g4hit) {

  PHG4Particle* particle = _truthinfo->GetHit( g4hit->get_trkid() );
  if (particle->get_primary_id() != (int)(0xFFFFFFFF)) {
    particle = _truthinfo->GetHit( particle->get_primary_id() );
  }
  
  return particle;
}

int CaloTruthEval::get_embed(PHG4Particle* particle) {
  
  return _truthinfo->isEmbeded(particle->get_track_id());
}

PHG4VtxPoint* CaloTruthEval::get_vertex(PHG4Particle* particle) {

  return _truthinfo->GetVtx( particle->get_vtx_id() );
}

bool CaloTruthEval::is_primary(PHG4Particle* particle) {

  if ((_do_cache) &&
      (_cache_is_primary.find(particle) != _cache_is_primary.end())) {
    return _cache_is_primary[particle];
  }
  
  bool is_primary = false;  
  PHG4TruthInfoContainer::Map primary_map = _truthinfo->GetPrimaryMap();
  for (PHG4TruthInfoContainer::ConstIterator iter = primary_map.begin(); 
       iter != primary_map.end(); 
       ++iter) {
    if (iter->second->get_track_id() == particle->get_track_id() ) {
      is_primary = true;
    }
  }

  if (_do_cache) _cache_is_primary.insert(make_pair(particle,is_primary));
  
  return is_primary;
}

std::set<PHG4Hit*> CaloTruthEval::get_shower_from_primary(PHG4Particle* primary) {

  if (!is_primary(primary)) return std::set<PHG4Hit*>();

  if ((_do_cache) &&
      (_cache_get_shower_from_primary.find(primary) != _cache_get_shower_from_primary.end())) {
    return _cache_get_shower_from_primary[primary];
  }
  
  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits
  for (PHG4HitContainer::ConstIterator g4iter = _g4hits->getHits().first;
       g4iter != _g4hits->getHits().second;
       ++g4iter) {

    PHG4Hit* g4hit = g4iter->second;
    PHG4Particle* candidate = get_primary_particle(g4hit);
    if (candidate->get_track_id() != primary->get_track_id()) continue;
    truth_hits.insert(g4hit);
  }

  if (_do_cache) _cache_get_shower_from_primary.insert(make_pair(primary,truth_hits));
  
  return truth_hits;
}

// moliere (90% containment) radius of scintilator hits
float CaloTruthEval::get_shower_moliere_radius(PHG4Particle* primary) {

  if (!is_primary(primary)) return NAN;

  if ((_do_cache) &&
      (_cache_get_shower_moliere_radius.find(primary) != _cache_get_shower_moliere_radius.end())) {
    return _cache_get_shower_moliere_radius[primary];
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

  if (shower_e == 0.0) return NAN; // no energy deposit!
  
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

  if (!is_primary(primary)) return NAN;

  if ((_do_cache) &&
      (_cache_get_shower_energy_deposit.find(primary) != _cache_get_shower_energy_deposit.end())) {
    return _cache_get_shower_energy_deposit[primary];
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

void CaloTruthEval::get_node_pointers(PHCompositeNode *topNode) {

  // need things off of the DST...
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }
  
  std::string name = "G4HIT_" + _caloname;
  _g4hits = findNode::getClass<PHG4HitContainer>(topNode,name.c_str());
  if (!_g4hits) {
    cerr << PHWHERE << " ERROR: Can't find " << name << endl;
    exit(-1);
  }
  
  return;
}
