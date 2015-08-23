
#ifndef __SVTXVERTEXEVAL_H__
#define __SVTXVERTEXEVAL_H__

#include "SvtxTrackEval.h"
#include "SvtxClusterEval.h"
#include "SvtxHitEval.h"

#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <set>
#include <map>

class SvtxVertexEval {

public:

  SvtxVertexEval(PHCompositeNode *topNode);
  virtual ~SvtxVertexEval() {}

  void next_event(PHCompositeNode *topNode);
  
  // access the sub evals (and the cached values)
  SvtxTrackEval*   get_track_eval() {return &_trackeval;}
  SvtxClusterEval* get_cluster_eval() {return _trackeval.get_cluster_eval();}
  SvtxHitEval*     get_hit_eval() {return _trackeval.get_hit_eval();}
  
  // backtrace through to PHG4Hits
  std::set<PHG4Particle*>  all_truth_particles (SvtxVertex* vertex);
  
  // backtrace through to PHG4VtxPoints
  std::set<PHG4VtxPoint*> all_truth_points           (SvtxVertex* vertex);
  PHG4VtxPoint*           max_truth_point_by_ntracks (SvtxVertex* vertex);

  // forwardtrace through to SvtxVertexs
  std::set<SvtxVertex*> all_vertexes_from (PHG4VtxPoint* truthpoint);
  SvtxVertex*           best_vertex_from  (PHG4VtxPoint* truthpoint);
  
  // overlap calculations
  unsigned int get_ntracks_contribution(SvtxVertex* svtxvertex, PHG4VtxPoint* truthpoint);  
  
private:

  void get_node_pointers(PHCompositeNode* topNode);
  
  SvtxTrackEval _trackeval;
  SvtxVertexMap* _vertexmap;
  SvtxTrackMap* _trackmap;
  PHG4TruthInfoContainer* _truthinfo;

  bool _do_cache;
  std::map<SvtxVertex*,std::set<PHG4Particle*> >               _cache_all_truth_particles;
  std::map<SvtxVertex*,std::set<PHG4VtxPoint*> >               _cache_all_truth_points;
  std::map<SvtxVertex*,PHG4VtxPoint*>                          _cache_max_truth_point_by_ntracks;
  std::map<PHG4VtxPoint*,std::set<SvtxVertex*> >               _cache_all_vertexes_from_point;
  std::map<PHG4VtxPoint*,SvtxVertex* >                         _cache_best_vertex_from_point;
  std::map<std::pair<SvtxVertex*,PHG4VtxPoint*>, unsigned int> _cache_get_ntracks_contribution;
};

#endif // __SVTXVERTEXEVAL_H__
