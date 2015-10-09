
#include "SvtxTruthEval.h"

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
#include <cassert>
#include <iostream>

using namespace std;

SvtxTruthEval::SvtxTruthEval(PHCompositeNode* topNode)
  : _truthinfo(NULL),
    _g4hits_svtx(NULL),
    _g4hits_tracker(NULL),
    _strict(true),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_hits_g4particle(),
    _cache_get_innermost_truth_hit(),
    _cache_get_outermost_truth_hit(),
    _cache_get_primary_g4hit() {
  get_node_pointers(topNode);
}

void SvtxTruthEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_hits_g4particle.clear();
  _cache_get_innermost_truth_hit.clear();
  _cache_get_outermost_truth_hit.clear();
  _cache_get_primary_g4hit.clear();
  
  get_node_pointers(topNode);
}

/// \todo this copy may be too expensive to call a lot...
std::set<PHG4Hit*> SvtxTruthEval::all_truth_hits() {

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

  if (_do_cache) _cache_all_truth_hits = truth_hits;

  return truth_hits;
}

std::set<PHG4Hit*> SvtxTruthEval::all_truth_hits(PHG4Particle* particle) {

  if (_strict) assert(particle);
  else if (!particle) return std::set<PHG4Hit*>();
  
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
  
  if (_do_cache) _cache_all_truth_hits_g4particle.insert(make_pair(particle,truth_hits));
  
  return truth_hits;
}

PHG4Hit* SvtxTruthEval::get_innermost_truth_hit(PHG4Particle* particle) {

  if (_strict) assert(particle);
  else if (!particle) return NULL;
  
  PHG4Hit* innermost_hit = NULL;
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

  if (_strict) assert(particle);
  else if (!particle) return NULL;
  
  PHG4Hit* outermost_hit = NULL;
  float outermost_radius = FLT_MAX*-1.0;
  
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

  return outermost_hit;
}

PHG4Particle* SvtxTruthEval::get_particle(PHG4Hit* g4hit) {

  if (_strict) assert(g4hit);
  else if (!g4hit) return NULL;
  
  PHG4Particle* particle = _truthinfo->GetHit( g4hit->get_trkid() );
  if (_strict) assert(particle); 

  return particle;
}

int SvtxTruthEval::get_embed(PHG4Particle* particle) {

  if (_strict) assert(particle);
  else if (!particle) return 0;

  return _truthinfo->isEmbeded(particle->get_track_id());
}

PHG4VtxPoint* SvtxTruthEval::get_vertex(PHG4Particle* particle) {

  if (_strict) assert(particle);
  else if (!particle) return NULL;

  if (particle->get_primary_id() == -1) {
    return _truthinfo->GetPrimaryVtx( particle->get_vtx_id() );  
  }

  PHG4VtxPoint* vtx = _truthinfo->GetVtx( particle->get_vtx_id() );
  if (_strict) assert(vtx);
 
  return vtx;
}

bool SvtxTruthEval::is_primary(PHG4Particle* particle) {

  if (_strict) assert(particle);
  else if (!particle) return false;
  
  bool is_primary = false;
  if (particle->get_parent_id() == 0) {
    is_primary = true;
  }
  
  return is_primary;
}

PHG4Particle* SvtxTruthEval::get_primary(PHG4Hit* g4hit) {

  if (_strict) assert(g4hit);
  else if (!g4hit) return NULL;

  if (_do_cache) {
    std::map<PHG4Hit*,PHG4Particle*>::iterator iter =
      _cache_get_primary_g4hit.find(g4hit);
    if (iter != _cache_get_primary_g4hit.end()) {
      return iter->second;
    }
  }
  
  PHG4Particle* particle = get_particle(g4hit);
  PHG4Particle* primary = get_primary(particle);

  if (_do_cache) _cache_get_primary_g4hit.insert(make_pair(g4hit,primary));
  
  if (_strict) assert(primary);
  
  return primary;
}

PHG4Particle* SvtxTruthEval::get_primary(PHG4Particle* particle) {

  if (_strict) assert(particle);
  else if (!particle) return NULL;

  // always report the primary from the Primary Map regardless if a
  // primary from the full Map was the argument
  PHG4Particle* primary = NULL;
  if (particle->get_primary_id() == -1) {
    primary = _truthinfo->GetPrimaryHit( particle->get_track_id() );
  } else {
    primary = _truthinfo->GetPrimaryHit( particle->get_primary_id() );
  }

  if (_strict) assert(primary);
  
  return primary;
}

 bool SvtxTruthEval::is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle) {

  if (_strict) {
    assert(g4hit);
    assert(particle);
  } else if (!g4hit||!particle) {
    return false;
  }

  bool is_from_particle = false;

  PHG4Particle* candidate = get_particle(g4hit);
  if (_strict) assert(candidate);

  if (candidate) {
    if (are_same_particle(candidate,particle)) {
      is_from_particle = true;
    }
  }

  return is_from_particle;
}

bool SvtxTruthEval::are_same_particle(PHG4Particle* p1, PHG4Particle* p2) {

  if (_strict) {
    assert(p1);
    assert(p2);
  } else if (!p1||!p2) {
    return false;
  }

  bool are_same_particle = false;

  if ((p1->get_primary_id() == -1) &&
      (p2->get_primary_id() == -1)) {
    // both particles are from the primary truth map, use track ids
    if (p1->get_track_id() == p2->get_track_id()) {
      are_same_particle = true;
    }
  } else if ((p1->get_primary_id() != -1) &&
	     (p2->get_primary_id() != -1)) {
    // both particles are from the total truth map, use track ids
    if (p1->get_track_id() == p2->get_track_id()) {
      are_same_particle = true;
    }   
  } else if ((p1->get_primary_id() == -1) &&
	     (p2->get_primary_id() != -1)) {
    // only the first particle is from the primary truth map, use one track id and one primary id
    if (p1->get_track_id() == p2->get_primary_id()) {
      are_same_particle = true;
    }
  } else if ((p1->get_primary_id() != -1) &&
	     (p2->get_primary_id() == -1)) {
    // only the second particle is from the primary truth map, use one track id and one primary id
    if (p1->get_primary_id() == p2->get_track_id()) {
      are_same_particle = true;
    }
  }
 
  return are_same_particle;
}
  
void SvtxTruthEval::get_node_pointers(PHCompositeNode* topNode) {

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  _g4hits_svtx    = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SVTX");
  _g4hits_tracker = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_SILICON_TRACKER");
  if (!_g4hits_svtx && !_g4hits_tracker) {
    cerr << PHWHERE << " ERROR: Can't find G4HIT_SVTX or G4HIT_SILICON_TRACKER" << endl;
    exit(-1);
  }
  
  return;
}
