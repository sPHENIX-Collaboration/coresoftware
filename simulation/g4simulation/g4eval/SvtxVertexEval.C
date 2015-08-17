
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
  : _topNode(topNode),
    _trackeval(topNode),
    _do_cache(true),
    _cache_all_truth_particles(),
    _cache_all_truth_points(),
    _cache_max_truth_point_by_ntracks(),
    _cache_all_vertexes_from_point(),
    _cache_best_vertex_from_point(),
    _cache_get_ntracks_contribution() {
}


void SvtxVertexEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_particles.clear();
  _cache_all_truth_points.clear();
  _cache_max_truth_point_by_ntracks.clear();
  _cache_all_vertexes_from_point.clear();
  _cache_best_vertex_from_point.clear();
  _cache_get_ntracks_contribution.clear();
  
  _trackeval.next_event(topNode);
  
  _topNode = topNode;  
}

std::set<PHG4Particle*> SvtxVertexEval::all_truth_particles(SvtxVertex* vertex) {

 if ((_do_cache) && (_cache_all_truth_particles.find(vertex) !=
		     _cache_all_truth_particles.end())) {
    return _cache_all_truth_particles[vertex];
  }
  
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(_topNode,"SvtxTrackMap");
  if (!trackmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
    exit(-1);
  }

  std::set<PHG4Particle*> all_particles;
  
  // loop over all tracks on vertex
  for (SvtxVertex::TrackIter iter = vertex->begin_tracks();
       iter != vertex->end_tracks();
       ++iter) {
    
    SvtxTrack* track = trackmap->get(*iter);
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
  
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(_topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }

  std::set<PHG4VtxPoint*> points;
  
  std::set<PHG4Particle*> particles = all_truth_particles(vertex);
  for (std::set<PHG4Particle*>::iterator iter = particles.begin();
       iter != particles.end();
       ++iter) {
    PHG4Particle* particle = *iter;
    points.insert(truthinfo->GetPrimaryVtx(particle->get_vtx_id()));
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
  
  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(_topNode,"SvtxVertexMap");
  if (!vertexmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap" << endl;
    exit(-1);
  }

  std::set<SvtxVertex*> all_vertexes;

  // loop over all vertexes on node
  for (SvtxVertexMap::Iter iter = vertexmap->begin();
       iter != vertexmap->end();
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
  
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(_topNode,"SvtxTrackMap");
  if (!trackmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
    exit(-1);
  }

  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(_topNode,"G4TruthInfo");
  if (!truthinfo) {
    cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
    exit(-1);
  }
  
  unsigned int ntracks = 0;

  for (SvtxVertex::TrackIter iter = vertex->begin_tracks();
       iter != vertex->end_tracks();
       ++iter) {
    
    SvtxTrack* track = trackmap->get(*iter);
    PHG4Particle* particle = _trackeval.max_truth_particle_by_nclusters(track);
    PHG4VtxPoint* candidate = truthinfo->GetPrimaryVtx(particle->get_vtx_id());
    if (candidate->get_id() == truthpoint->get_id()) ++ntracks;
  }
  
  if (_do_cache) _cache_get_ntracks_contribution.insert(make_pair(make_pair(vertex,truthpoint),ntracks));
  
  return ntracks;
}
