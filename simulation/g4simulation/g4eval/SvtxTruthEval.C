#include "SvtxTruthEval.h"

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
#include <cassert>
#include <iostream>

using namespace std;

SvtxTruthEval::SvtxTruthEval(PHCompositeNode* topNode)
  : _basetrutheval(topNode),
    _truthinfo(nullptr),
    _g4hits_svtx(nullptr),
    _g4hits_tracker(nullptr),
    _g4hits_maps(nullptr),
    _strict(false),
    _verbosity(1),
    _errors(0), 
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_hits_g4particle(),
    _cache_get_innermost_truth_hit(),
    _cache_get_outermost_truth_hit(),
    _cache_get_primary_particle_g4hit() {
  get_node_pointers(topNode);
}

SvtxTruthEval::~SvtxTruthEval() {
  if (_verbosity > 0) {
    if ((_errors > 0)||(_verbosity > 1)) {
      cout << "SvtxTruthEval::~SvtxTruthEval() - Error Count: " << _errors << endl;
    }
  }
}

void SvtxTruthEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_hits_g4particle.clear();
  _cache_get_innermost_truth_hit.clear();
  _cache_get_outermost_truth_hit.clear();
  _cache_get_primary_particle_g4hit.clear();

  _basetrutheval.next_event(topNode);
  
  get_node_pointers(topNode);
}

/// \todo this copy may be too expensive to call a lot...
std::set<PHG4Hit*> SvtxTruthEval::all_truth_hits() {

  if (!has_node_pointers()) {++_errors; return std::set<PHG4Hit*>();}
  
  if (_do_cache) {
    if (!_cache_all_truth_hits.empty()) {
      return _cache_all_truth_hits;
    }
  }
  
  // since the SVTX can be composed of two different trackers this is a
  // handy function to spill out all the g4hits from both "detectors"
  
  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits in the cylinder layers
  if (_g4hits_svtx) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_svtx->getHits().first;
	 g4iter != _g4hits_svtx->getHits().second;
	 ++g4iter) {

      PHG4Hit* g4hit = g4iter->second;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the ladder layers
  if (_g4hits_tracker) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_tracker->getHits().first;
	 g4iter != _g4hits_tracker->getHits().second;
	 ++g4iter) {
      
      PHG4Hit* g4hit = g4iter->second;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the maps layers
  if (_g4hits_maps) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_maps->getHits().first;
	 g4iter != _g4hits_maps->getHits().second;
	 ++g4iter) {
      
      PHG4Hit* g4hit = g4iter->second;
      truth_hits.insert(g4hit);
    }
  }

  if (_do_cache) _cache_all_truth_hits = truth_hits;

  return truth_hits;
}

std::set<PHG4Hit*> SvtxTruthEval::all_truth_hits(PHG4Particle* particle) {

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

  // loop over all the g4hits in the cylinder layers
  if (_g4hits_svtx) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_svtx->getHits().first;
	 g4iter != _g4hits_svtx->getHits().second;
	 ++g4iter) {

      PHG4Hit* g4hit = g4iter->second;
      if (!is_g4hit_from_particle(g4hit,particle)) continue;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the ladder layers
  if (_g4hits_tracker) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_tracker->getHits().first;
	 g4iter != _g4hits_tracker->getHits().second;
	 ++g4iter) {
      
      PHG4Hit* g4hit = g4iter->second;
      if (!is_g4hit_from_particle(g4hit,particle)) continue;
      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the maps ladder layers
  if (_g4hits_maps) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_maps->getHits().first;
	 g4iter != _g4hits_maps->getHits().second;
	 ++g4iter) {
      
      PHG4Hit* g4hit = g4iter->second;
      if (!is_g4hit_from_particle(g4hit,particle)) continue;
      truth_hits.insert(g4hit);
    }
  }
  
  if (_do_cache) _cache_all_truth_hits_g4particle.insert(make_pair(particle,truth_hits));
  
  return truth_hits;
}

PHG4Hit* SvtxTruthEval::get_innermost_truth_hit(PHG4Particle* particle) {

  if (!has_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return nullptr;}
  
  PHG4Hit* innermost_hit = nullptr;
  float innermost_radius = FLT_MAX;
  
  std::set<PHG4Hit*> truth_hits = all_truth_hits(particle);
  for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
       iter != truth_hits.end();
       ++iter) {
    PHG4Hit* candidate = *iter;
    float x = candidate->get_x(0); // use entry points
    float y = candidate->get_y(0); // use entry points
    float r = sqrt(x*x+y*y);
    if (r < innermost_radius) {
      innermost_radius = r;
      innermost_hit = candidate;      
    }
  }

  return innermost_hit;
}

PHG4Hit* SvtxTruthEval::get_outermost_truth_hit(PHG4Particle* particle) {

  if (!has_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(particle);}
  else if (!particle) {++_errors; return nullptr;}
  
  PHG4Hit* outermost_hit = nullptr;
  float outermost_radius = FLT_MAX*-1.0;
  
  if (_do_cache) {
    std::map<PHG4Particle*,PHG4Hit*>::iterator iter =
      _cache_get_outermost_truth_hit.find(particle);
    if (iter != _cache_get_outermost_truth_hit.end()) {
      return iter->second;
    }
  }

  std::set<PHG4Hit*> truth_hits = all_truth_hits(particle);
  for (std::set<PHG4Hit*>::iterator iter = truth_hits.begin();
       iter != truth_hits.end();
       ++iter) {
    PHG4Hit* candidate = *iter;
    float x = candidate->get_x(1); // use exit points
    float y = candidate->get_y(1); // use exit points
    float r = sqrt(x*x+y*y);
    if (r > outermost_radius) {
      outermost_radius = r;
      outermost_hit = candidate;      
    }
  }
  if (_do_cache) _cache_get_outermost_truth_hit.insert(make_pair(particle,outermost_hit));

  return outermost_hit;
}

PHG4Particle* SvtxTruthEval::get_particle(PHG4Hit* g4hit) {
  return _basetrutheval.get_particle(g4hit);
}

int SvtxTruthEval::get_embed(PHG4Particle* particle) {
  return _basetrutheval.get_embed(particle);
}

PHG4VtxPoint* SvtxTruthEval::get_vertex(PHG4Particle* particle) {
  return _basetrutheval.get_vertex(particle);
}

bool SvtxTruthEval::is_primary(PHG4Particle* particle) {
  return _basetrutheval.is_primary(particle);
}

PHG4Particle* SvtxTruthEval::get_primary_particle(PHG4Hit* g4hit) {

  if (!has_node_pointers()) {++_errors; return nullptr;}
  
  if (_strict) {assert(g4hit);}
  else if (!g4hit) {++_errors; return nullptr;}

  if (_do_cache) {
    std::map<PHG4Hit*,PHG4Particle*>::iterator iter =
      _cache_get_primary_particle_g4hit.find(g4hit);
    if (iter != _cache_get_primary_particle_g4hit.end()) {
      return iter->second;
    }
  }
  
  PHG4Particle* primary = _basetrutheval.get_primary_particle(g4hit);

  if (_do_cache) _cache_get_primary_particle_g4hit.insert(make_pair(g4hit,primary));
  
  if (_strict) {assert(primary);}
  else if (!primary) {++_errors;}
  
  return primary;
}

PHG4Particle* SvtxTruthEval::get_primary_particle(PHG4Particle* particle) {
  return _basetrutheval.get_primary_particle(particle);
}

bool SvtxTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle) {
  return _basetrutheval.is_g4hit_from_particle(g4hit,particle);
}

bool SvtxTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2) {
  return _basetrutheval.are_same_particle(p1,p2);
}

bool SvtxTruthEval::are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2) {
  return _basetrutheval.are_same_vertex(vtx1,vtx2);
}
  
void SvtxTruthEval::get_node_pointers(PHCompositeNode* topNode) {
  
  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");

  _g4hits_svtx    = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_TPC");
  _g4hits_tracker = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SILICON_TRACKER");
  _g4hits_maps = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_MAPS");

  return;
}

bool SvtxTruthEval::has_node_pointers() {

  if (_strict) assert(_truthinfo);
  else if (!_truthinfo) return false;

  if (_strict) assert(_g4hits_svtx || _g4hits_tracker || _g4hits_maps);
  else if (!_g4hits_svtx && !_g4hits_tracker && !_g4hits_maps) return false;
  
  return true;
}
