#include "SvtxHitEval.h"

#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <set>
#include <map>
#include <float.h>
#include <cassert>

using namespace std;

SvtxHitEval::SvtxHitEval(PHCompositeNode* topNode)
  : _trutheval(topNode),
    _hitmap(nullptr),
    _g4cells_svtx(nullptr),
    _g4cells_tracker(nullptr),
    _g4cells_maps(nullptr),
    _g4hits_svtx(nullptr),
    _g4hits_tracker(nullptr),
   _g4hits_maps(nullptr),
    _truthinfo(nullptr),
    _strict(false),
    _verbosity(1),
    _errors(0), 
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_max_truth_hit_by_energy(),
    _cache_all_truth_particles(),
    _cache_max_truth_particle_by_energy(),
    _cache_all_hits_from_particle(),
    _cache_all_hits_from_g4hit(),
    _cache_best_hit_from_g4hit(),
    _cache_get_energy_contribution_g4particle(),
    _cache_get_energy_contribution_g4hit() {
  get_node_pointers(topNode);
}

SvtxHitEval::~SvtxHitEval() {
  if (_verbosity > 0) {
    if ((_errors > 0)||(_verbosity > 1)) {
      cout << "SvtxHitEval::~SvtxHitEval() - Error Count: " << _errors << endl;
    }
  }
}

void SvtxHitEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_max_truth_hit_by_energy.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_energy.clear();
  _cache_all_hits_from_particle.clear();
  _cache_all_hits_from_g4hit.clear();
  _cache_best_hit_from_g4hit.clear();
  _cache_get_energy_contribution_g4particle.clear();
  _cache_get_energy_contribution_g4hit.clear();

  _trutheval.next_event(topNode);
  
  get_node_pointers(topNode);
}

PHG4Cell* SvtxHitEval::get_cell(SvtxHit* hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  if (_strict) {assert(hit);}
  else if (!hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  // hop from reco hit to g4cell
  PHG4Cell* cell = nullptr;
  if (!cell&&_g4cells_svtx)    cell = _g4cells_svtx->findCell(hit->get_cellid());
  if (!cell&&_g4cells_tracker) cell = _g4cells_tracker->findCell(hit->get_cellid());
  if (!cell&&_g4cells_maps) cell = _g4cells_maps->findCell(hit->get_cellid());

  // only noise hits (cellid left at default value) should not trace
  if ((_strict) && (hit->get_cellid() != 0xFFFFFFFF)) {assert(cell);}
  else if (!cell) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl;}
  
  return cell;
}

std::set<PHG4Hit*> SvtxHitEval::all_truth_hits(SvtxHit* hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<PHG4Hit*>();}
  
  if (_strict) {assert(hit);}
  else if (!hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<PHG4Hit*>();}
  
  if (_do_cache) {
    std::map<SvtxHit*,std::set<PHG4Hit*> >::iterator iter =
      _cache_all_truth_hits.find(hit);
    if (iter != _cache_all_truth_hits.end()) {
      return iter->second;
    }
  }
    
  std::set<PHG4Hit*> truth_hits;
  
  // hop from reco hit to g4cell
  PHG4Cell *cell = nullptr;
  if (!cell&&_g4cells_svtx)    cell = _g4cells_svtx->findCell(hit->get_cellid());
  if (!cell&&_g4cells_tracker) cell = _g4cells_tracker->findCell(hit->get_cellid());
  if (!cell&&_g4cells_maps) cell = _g4cells_maps->findCell(hit->get_cellid());

  // only noise hits (cellid left at default value) should not trace
  if ((_strict) && (hit->get_cellid() != 0xFFFFFFFF)) assert(cell);
  else if (!cell) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return truth_hits;}

  //cout << "Eval: hitid " << hit->get_id() << " cellid " << cell->get_cellid() << endl;
  // loop over all the g4hits in this cell
  for (PHG4Cell::EdepConstIterator g4iter = cell->get_g4hits().first;
       g4iter != cell->get_g4hits().second;
       ++g4iter) {
    //cout << "    Looking for hit " << g4iter->first << " in layer " << cell->get_layer() << " with edep " << g4iter->second << endl;      
    PHG4Hit* g4hit = nullptr;
    if (!g4hit&&_g4hits_svtx)    g4hit = _g4hits_svtx->findHit(g4iter->first);
    if (!g4hit&&_g4hits_tracker) g4hit = _g4hits_tracker->findHit(g4iter->first);
    if (!g4hit&&_g4hits_maps) g4hit = _g4hits_maps->findHit(g4iter->first);
    if(!g4hit) cout << "    Failed to find  g4hit " << g4iter->first << " with edep " << g4iter->second << endl;
    if (_strict) assert(g4hit);
    else if (!g4hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; continue;}
    
    // fill output set
    truth_hits.insert(g4hit);
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(hit,truth_hits));
  
  return truth_hits;
}

PHG4Hit* SvtxHitEval::max_truth_hit_by_energy(SvtxHit* hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  if (_strict) {assert(hit);}
  else if (!hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  if (_do_cache) {
    std::map<SvtxHit*,PHG4Hit*>::iterator iter =
      _cache_max_truth_hit_by_energy.find(hit);
    if (iter != _cache_max_truth_hit_by_energy.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> hits = all_truth_hits(hit);
  PHG4Hit* max_hit = nullptr;
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

  if (_do_cache) _cache_max_truth_hit_by_energy.insert(make_pair(hit,max_hit));
  
  return max_hit;
}
  
std::set<PHG4Particle*> SvtxHitEval::all_truth_particles(SvtxHit* hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<PHG4Particle*>();}
  
  if (_strict) {assert(hit);}
  else if (!hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<PHG4Particle*>();}
  
  if (_do_cache) {
    std::map<SvtxHit*,std::set<PHG4Particle*> >::iterator iter =
      _cache_all_truth_particles.find(hit);
    if (iter != _cache_all_truth_particles.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Particle*> truth_particles;
  
  std::set<PHG4Hit*> g4hits = all_truth_hits(hit);

  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* g4hit = *iter;
    PHG4Particle* particle = get_truth_eval()->get_particle( g4hit );

    if (_strict) assert(particle);
    else if (!particle) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; continue;}

    truth_particles.insert(particle);
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(hit,truth_particles));
  
  return truth_particles;
}

PHG4Particle* SvtxHitEval::max_truth_particle_by_energy(SvtxHit* hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  if (_strict) {assert(hit);}
  else if (!hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  if (_do_cache) {
    std::map<SvtxHit*,PHG4Particle*>::iterator iter =
      _cache_max_truth_particle_by_energy.find(hit);
    if (iter != _cache_max_truth_particle_by_energy.end()) {
      return iter->second;
    }
  }
  
  // loop over all particles associated with this hit and
  // get the energy contribution for each one, record the max
  PHG4Particle* max_particle = nullptr;
  float max_e = FLT_MAX*-1.0;
  std::set<PHG4Particle*> particles = all_truth_particles(hit);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter) {

    PHG4Particle* particle = *iter;
    float e = get_energy_contribution(hit,particle);
    if (e > max_e) {
      max_e = e;
      max_particle = particle;      
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_energy.insert(make_pair(hit,max_particle));
  
  return max_particle;
}

std::set<SvtxHit*> SvtxHitEval::all_hits_from(PHG4Particle* g4particle) { 

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<SvtxHit*>();}
  
  if (_strict) {assert(g4particle);}
  else if (!g4particle) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<SvtxHit*>();}
  
  if (_do_cache) {
    std::map<PHG4Particle*,std::set<SvtxHit*> >::iterator iter =
      _cache_all_hits_from_particle.find(g4particle);
    if (iter != _cache_all_hits_from_particle.end()) {
      return iter->second;
    }
  }
 
  std::set<SvtxHit*> hits;
  
  // loop over all the hits
  for (SvtxHitMap::Iter iter = _hitmap->begin();
       iter != _hitmap->end();
       ++iter) {

    SvtxHit* hit = iter->second;
    // loop over all truth particles connected to this hit
    std::set<PHG4Particle*> g4particles = all_truth_particles(hit);
    for (std::set<PHG4Particle*>::iterator jter = g4particles.begin();
	 jter != g4particles.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;
      if (get_truth_eval()->are_same_particle(candidate,g4particle)) {
	hits.insert(hit);
      }    
    }
  }

  if (_do_cache) _cache_all_hits_from_particle.insert(make_pair(g4particle,hits));
  
  return hits;
}

std::set<SvtxHit*> SvtxHitEval::all_hits_from(PHG4Hit* g4hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<SvtxHit*>();}
  
  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return std::set<SvtxHit*>();}
  
  if (_do_cache) {
    std::map<PHG4Hit*,std::set<SvtxHit*> >::iterator iter =
      _cache_all_hits_from_g4hit.find(g4hit);
    if (iter != _cache_all_hits_from_g4hit.end()) {
      return iter->second;
    }
  }
  
  std::set<SvtxHit*> hits;

   unsigned int hit_layer = g4hit->get_layer();
  
  // loop over all the hits
  for (SvtxHitMap::Iter iter = _hitmap->begin();
       iter != _hitmap->end();
       ++iter) {

    SvtxHit* hit = iter->second;

    if (hit->get_layer() != hit_layer) continue;
    
    // loop over all truth hits connected to this hit
    std::set<PHG4Hit*> g4hits = all_truth_hits(hit);
    for (std::set<PHG4Hit*>::iterator jter = g4hits.begin();
	 jter != g4hits.end();
	 ++jter) {
      PHG4Hit* candidate = *jter;
      if (candidate->get_hit_id() == g4hit->get_hit_id()) {
	hits.insert(hit);
      }    
    }
  }

  if (_do_cache) _cache_all_hits_from_g4hit.insert(make_pair(g4hit,hits));
  
  return hits;
}

SvtxHit* SvtxHitEval::best_hit_from(PHG4Hit* g4hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return nullptr;}
  
  if (_do_cache) {
    std::map<PHG4Hit*,SvtxHit*>::iterator iter =
      _cache_best_hit_from_g4hit.find(g4hit);
    if (iter != _cache_best_hit_from_g4hit.end()) {
      return iter->second;
    }
  }

  SvtxHit* best_hit = nullptr;
  float best_energy = 0.0;  
  std::set<SvtxHit*> hits = all_hits_from(g4hit);
  for (std::set<SvtxHit*>::iterator iter = hits.begin();
       iter != hits.end();
       ++iter) {
    SvtxHit* hit = *iter;
    float energy = get_energy_contribution(hit,g4hit);
    if (energy > best_energy) {
      best_hit = hit;
      best_energy = energy;
    }
  }
 
  if (_do_cache) _cache_best_hit_from_g4hit.insert(make_pair(g4hit,best_hit));
  
  return best_hit;
}

// overlap calculations
float SvtxHitEval::get_energy_contribution(SvtxHit* hit, PHG4Particle* particle) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return NAN;}
  
  if (_strict) {
    assert(hit);
    assert(particle);
  } else if (!hit||!particle) {
    ++_errors; cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }
  
  if (_do_cache) {
    std::map<std::pair<SvtxHit*,PHG4Particle*>,float>::iterator iter =
      _cache_get_energy_contribution_g4particle.find(make_pair(hit,particle));
    if (iter != _cache_get_energy_contribution_g4particle.end()) {
      return iter->second;
    }
  }
  
  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hit);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* g4hit = *iter;
    if (get_truth_eval()->is_g4hit_from_particle(g4hit,particle)) {
      energy += g4hit->get_edep();
    }
  }

  if (_do_cache) _cache_get_energy_contribution_g4particle.insert(make_pair(make_pair(hit,particle),energy));
  
  return energy;
}

float SvtxHitEval::get_energy_contribution(SvtxHit* hit, PHG4Hit* g4hit) {

  if (!has_node_pointers()) {++_errors; cout << PHWHERE << " nerr: " << _errors << endl; return NAN;}
  
  if (_strict) {
    assert(hit);
    assert(g4hit);
  } else if (!hit||!g4hit) {
    ++_errors; cout << PHWHERE << " nerr: " << _errors << endl;
    return NAN;
  }
  
  if (_do_cache) {
    std::map<std::pair<SvtxHit*,PHG4Hit*>,float>::iterator iter =
      _cache_get_energy_contribution_g4hit.find(make_pair(hit,g4hit));
    if (iter != _cache_get_energy_contribution_g4hit.end()) {
      return iter->second;
    }
  }
    
  // this is a fairly simple existance check right now, but might be more
  // complex in the future, so this is here mostly as future-proofing.
  
  float energy = 0.0;
  std::set<PHG4Hit*> g4hits = all_truth_hits(hit);
  for (std::set<PHG4Hit*>::iterator iter = g4hits.begin();
       iter != g4hits.end();
       ++iter) {
    PHG4Hit* candidate = *iter;
    if (candidate->get_hit_id() != g4hit->get_hit_id()) continue;  
    energy += candidate->get_edep();
  }

  if (_do_cache) _cache_get_energy_contribution_g4hit.insert(make_pair(make_pair(hit,g4hit),energy));
  
  return energy;
}

void SvtxHitEval::get_node_pointers(PHCompositeNode* topNode) {

  // need things off of the DST...
  _hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");

  // need things off of the DST...
  _g4cells_svtx    = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SVTX");
  _g4cells_tracker = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SILICON_TRACKER");
  _g4cells_maps = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_MAPS");

  //  _g4hits_svtx    = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SVTX");
  _g4hits_svtx    = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_TPC");
  _g4hits_tracker = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SILICON_TRACKER");
  _g4hits_maps = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_MAPS");
  
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  
  return;
}

bool SvtxHitEval::has_node_pointers() {

  if (_strict) assert(_hitmap);
  else if (!_hitmap) 
    {
      return false;
    }

  if (_strict) assert(_g4cells_svtx || _g4cells_tracker || _g4cells_maps);
  else if (!_g4cells_svtx && !_g4cells_tracker && !_g4cells_maps)
    { 
      cout << "no cells" << endl;
      return false;
    }

  if (_strict) assert(_g4hits_svtx || _g4hits_tracker || _g4hits_maps);
  else if (!_g4hits_svtx && !_g4hits_tracker && !_g4hits_maps) 
    {
      cout << "no hits" << endl;
      return false; 
    }
  if (_strict) assert(_truthinfo);
  else if (!_truthinfo)
    {
      cout << " no truth" << endl;
      return false;
    }

  return true;
}
