
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

using namespace std;

SvtxTruthEval::SvtxTruthEval(PHCompositeNode* topNode)
  : _truthinfo(NULL),
    _g4hits_svtx(NULL),
    _g4hits_tracker(NULL),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_hits_g4particle(),
    _cache_get_innermost_truth_hit(),
    _cache_get_outermost_truth_hit() {
  get_node_pointers(topNode);
}

void SvtxTruthEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_hits_g4particle.clear();
  _cache_get_innermost_truth_hit.clear();
  _cache_get_outermost_truth_hit.clear();
  
  get_node_pointers(topNode);
}

/// \todo this copy may be too expensive to call a lot...
std::set<PHG4Hit*> SvtxTruthEval::all_truth_hits() {

  if ((_do_cache)&&(!_cache_all_truth_hits.empty())) return _cache_all_truth_hits;
  
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

  if ((_do_cache) &&
      (_cache_all_truth_hits_g4particle.find(particle) != _cache_all_truth_hits_g4particle.end())) {
    return _cache_all_truth_hits_g4particle[particle];
  }
  
  std::set<PHG4Hit*> truth_hits;

  // loop over all the g4hits in the cylinder layers
  if (_g4hits_svtx) {
    for (PHG4HitContainer::ConstIterator g4iter = _g4hits_svtx->getHits().first;
	 g4iter != _g4hits_svtx->getHits().second;
	 ++g4iter) {

      PHG4Hit* g4hit = g4iter->second;
      if (g4hit->get_trkid() != particle->get_track_id()) continue;
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

  PHG4Hit* outermost_hit = NULL;
  float outermost_radius = FLT_MAX-1.0;
  
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

  PHG4Particle* particle = _truthinfo->GetHit( g4hit->get_trkid() );
  return particle;
}

int SvtxTruthEval::get_embed(PHG4Particle* particle) {

  return _truthinfo->isEmbeded(particle->get_track_id());
}

bool SvtxTruthEval::is_primary(PHG4Particle* particle) {

  bool is_primary = false;  
  PHG4TruthInfoContainer::Map primary_map = _truthinfo->GetPrimaryMap();
  for (PHG4TruthInfoContainer::ConstIterator iter = primary_map.begin(); 
       iter != primary_map.end(); 
       ++iter) {
    if (iter->second->get_track_id() == particle->get_track_id() ) {
      is_primary = true;
    }
  }
  
  return is_primary;
}

PHG4VtxPoint* SvtxTruthEval::get_vertex(PHG4Particle* particle) {

  return _truthinfo->GetVtx( particle->get_vtx_id() );
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
