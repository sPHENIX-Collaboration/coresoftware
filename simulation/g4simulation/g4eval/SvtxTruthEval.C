
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

      // if particle is from the primary map it may have a different id
      // than the copy of the particle that left the hit
      // if (particle->get_primary_id() == -1) {
      // 	if (g4hit->get
      // }
      if (g4hit->get_trkid() != particle->get_track_id()) continue;

      // ********

      truth_hits.insert(g4hit);
    }
  }

  // loop over all the g4hits in the ladder layers
  if (_g4hits_tracker) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_tracker->getHits().first;
	 g4iter != _g4hits_tracker->getHits().second;
	 ++g4iter) {
      
      PHG4Hit* g4hit = g4iter->second;
      if (g4hit->get_trkid() != particle->get_track_id()) continue;
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

  // we have three cases to test:
  // #1 - input particle is not primary and resides in the truthinfo map
  //        in this case the track id on the hit must match directly for true return
  // #2 - input particle is primary and resides in the truthinfo map
  //        in this case the track id on the hit must match directly for true return
  // #3 - input particle is primary and resides in the truthinfo primary map
  //        in this case the track id on the hit will match to the copy id in the truthinfo map
  //        the particle leaving the g4hit must be both a primary (parent id==0) and
  //        must have a primary particle id that matches our input particle track id
  
  bool is_from_particle = false;
  if (particle->get_primary_id() != -1) {
    // particle came from the total truth map
    // and so will leave a hit with a matching id
    if (g4hit->get_trkid() == particle->get_track_id()) {
      is_from_particle = true;
    }
  } else {
    // particle came from the primary truth map
    // if it left a hit that will trace directly to the copy with a possibly different id
    // but the particle that left the g4hit will think it is a primary
    // and will trace to the primary particle map entry
    
    PHG4Particle* candidate = get_particle(g4hit);
    if (_strict) assert(candidate);

    if (candidate) {
      if (is_primary(candidate)) {
	if (particle->get_track_id() == candidate->get_primary_id()) {
	  is_from_particle = true;
	}
      }
    }
  }

  return is_from_particle;
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
