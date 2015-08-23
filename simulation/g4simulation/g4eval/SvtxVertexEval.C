
#include "SvtxVertexEval.h"

#include "SvtxTrackEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <cstdlib>
#include <set>
#include <float.h>
#include <algorithm>

using namespace std;

SvtxVertexEval::SvtxVertexEval(PHCompositeNode* topNode)
  : _trackeval(topNode),
    _vertexmap(NULL),
    _trackmap(NULL),
    _truthinfo(NULL),
    _do_cache(true),
    _cache_all_truth_particles(),
    _cache_all_truth_points(),
    _cache_max_truth_point_by_ntracks(),
    _cache_all_vertexes_from_point(),
    _cache_best_vertex_from_point(),
    _cache_get_ntracks_contribution() {
  get_node_pointers(topNode);
}


void SvtxVertexEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_particles.clear();
  _cache_all_truth_points.clear();
  _cache_max_truth_point_by_ntracks.clear();
  _cache_all_vertexes_from_point.clear();
  _cache_best_vertex_from_point.clear();
  _cache_get_ntracks_contribution.clear();
  
  _trackeval.next_event(topNode);
  
  get_node_pointers(topNode);
}

std::set<PHG4Particle*> SvtxVertexEval::all_truth_particles(SvtxVertex* vertex) {

 if ((_do_cache) && (_cache_all_truth_particles.find(vertex) !=
		     _cache_all_truth_particles.end())) {
    return _cache_all_truth_particles[vertex];
  }
 
  std::set<PHG4Particle*> all_particles;
  
  // loop over all tracks on vertex
  for (SvtxVertex::TrackIter iter = vertex->begin_tracks();
       iter != vertex->end_tracks();
       ++iter) {
    
    SvtxTrack* track = _trackmap->get(*iter);    
    all_particles.insert(_trackeval.max_truth_particle_by_nclusters(track));
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(vertex,all_particles));
  
  return all_particles;
}

std::set<PHG4VtxPoint*> SvtxVertexEval::all_truth_points(SvtxVertex* vertex) {

 if ((_do_cache) && (_cache_all_truth_points.find(vertex) !=
		     _cache_all_truth_points.end())) {
    return _cache_all_truth_points[vertex];
  }

  std::set<PHG4VtxPoint*> points;
  
  std::set<PHG4Particle*> particles = all_truth_particles(vertex);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter) {
    PHG4Particle* particle = *iter;

    // only consider primary particles
    PHG4TruthInfoContainer::Map primarymap = _truthinfo->GetPrimaryMap();
    if (primarymap.find(particle->get_track_id()) == primarymap.end()) continue;
        
    points.insert(_truthinfo->GetPrimaryVtx(particle->get_vtx_id()));
  }

  if (_do_cache) _cache_all_truth_points.insert(make_pair(vertex,points));
  
  return points;
}

PHG4VtxPoint* SvtxVertexEval::max_truth_point_by_ntracks(SvtxVertex* vertex) {

  if ((_do_cache) && (_cache_max_truth_point_by_ntracks.find(vertex) !=
		      _cache_max_truth_point_by_ntracks.end())) {
    return _cache_max_truth_point_by_ntracks[vertex];
  }
  
  std::set<PHG4VtxPoint*> points = all_truth_points(vertex);

  PHG4VtxPoint* max_point = NULL;
  unsigned int max_ntracks = 0;
  
  for (std::set<PHG4VtxPoint*>::iterator iter = points.begin();
      iter != points.end();
      ++iter) {
    PHG4VtxPoint* candidate = *iter;
    unsigned int ntracks = get_ntracks_contribution(vertex,candidate);
    if (ntracks > max_ntracks) {
      max_ntracks = ntracks;
      max_point = candidate;
    }
  }

  if (_do_cache) _cache_max_truth_point_by_ntracks.insert(make_pair(vertex,max_point));
  
  return max_point;
}
   
std::set<SvtxVertex*> SvtxVertexEval::all_vertexes_from(PHG4VtxPoint* truthpoint) {

 if ((_do_cache) && (_cache_all_vertexes_from_point.find(truthpoint) !=
		     _cache_all_vertexes_from_point.end())) {
    return _cache_all_vertexes_from_point[truthpoint];
  }

  std::set<SvtxVertex*> all_vertexes;

  // loop over all vertexes on node
  for (SvtxVertexMap::Iter iter = _vertexmap->begin();
       iter != _vertexmap->end();
       ++iter) {
    SvtxVertex* vertex = &iter->second;
    std::set<PHG4VtxPoint*> points = all_truth_points(vertex);
    for (std::set<PHG4VtxPoint*>::iterator jter = points.begin();
	 jter != points.end();
	 ++jter) {
      PHG4VtxPoint* point = *jter;
      if (point->get_id() == truthpoint->get_id()) {
	all_vertexes.insert(vertex);
      }
    }
  }
  
  if (_do_cache) _cache_all_vertexes_from_point.insert(make_pair(truthpoint,all_vertexes));
  
  return all_vertexes;
}

SvtxVertex* SvtxVertexEval::best_vertex_from(PHG4VtxPoint* truthpoint) {
 
  if ((_do_cache) && (_cache_best_vertex_from_point.find(truthpoint) !=
		      _cache_best_vertex_from_point.end())) {
    return _cache_best_vertex_from_point[truthpoint];
  }

  SvtxVertex* best_vertex = NULL;
  unsigned int best_count = 0;
  std::set<SvtxVertex*> tracks = all_vertexes_from(truthpoint);
  for (std::set<SvtxVertex*>::iterator iter = tracks.begin();
       iter != tracks.end();
       ++iter) {
    SvtxVertex* vertex = *iter;
    unsigned int count = get_ntracks_contribution(vertex,truthpoint);
    if (count > best_count) {
      best_vertex = vertex;
      best_count = count;
    }
  }
  
  if (_do_cache) _cache_best_vertex_from_point.insert(make_pair(truthpoint,best_vertex));
  
  return best_vertex;
}

// overlap calculations
unsigned int SvtxVertexEval::get_ntracks_contribution(SvtxVertex* vertex, PHG4VtxPoint* truthpoint) {

  if ((_do_cache) && (_cache_get_ntracks_contribution.find(make_pair(vertex,truthpoint)) !=
		      _cache_get_ntracks_contribution.end())) {
    return _cache_get_ntracks_contribution[make_pair(vertex,truthpoint)];
  }
  
  unsigned int ntracks = 0;

  for (SvtxVertex::TrackIter iter = vertex->begin_tracks();
       iter != vertex->end_tracks();
       ++iter) {
    
    SvtxTrack* track = _trackmap->get(*iter);
    PHG4Particle* particle = _trackeval.max_truth_particle_by_nclusters(track);

    // only consider primary particles
    PHG4TruthInfoContainer::Map primarymap = _truthinfo->GetPrimaryMap();
    if (primarymap.find(particle->get_track_id()) == primarymap.end()) continue;

    PHG4VtxPoint* candidate = _truthinfo->GetPrimaryVtx(particle->get_vtx_id());
    if (candidate->get_id() == truthpoint->get_id()) ++ntracks;
  }
  
  if (_do_cache) _cache_get_ntracks_contribution.insert(make_pair(make_pair(vertex,truthpoint),ntracks));
  
  return ntracks;
}

void SvtxVertexEval::get_node_pointers(PHCompositeNode* topNode) {

  // need things off the DST...
  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (!_vertexmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap" << endl;
    exit(-1);
  }
  
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if (!_trackmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
    exit(-1);
  }

  _truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }
  
  return;
}
