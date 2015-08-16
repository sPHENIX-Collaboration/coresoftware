
#include "SvtxVertexEval.h"

#include "SvtxTrackEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxVertexMap.h>
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
  return std::set<PHG4Particle*>();
}

std::set<PHG4VtxPoint*> SvtxVertexEval::all_truth_points(SvtxVertex* vertex) {
  return std::set<PHG4VtxPoint*>();
}

PHG4VtxPoint* SvtxVertexEval::max_truth_point_by_ntracks(SvtxVertex* vertex) {
  return NULL;
}
  
std::set<SvtxVertex*> SvtxVertexEval::all_vertexes_from(PHG4VtxPoint* truthpoint) {
  return std::set<SvtxVertex*>();
}

SvtxVertex* SvtxVertexEval::best_vertex_from(PHG4VtxPoint* truthpoint) {
  return NULL;
}

// overlap calculations
unsigned int SvtxVertexEval::get_ntracks_contribution(SvtxVertex* track, PHG4VtxPoint* truthpoint) {
  return 0;
}
