
#ifndef __SVTXEVALSTACK_H__
#define __SVTXEVALSTACK_H__

#include "SvtxTrackEval.h"
#include "SvtxClusterEval.h"
#include "SvtxHitEval.h"
#include "SvtxTruthEval.h"

#include <phool/PHCompositeNode.h>

// This user class provides pointers to the
// full set of tracking evaluators and
// protects the user from future introduction
// of new eval heirachies (new eval objects can
// be introduced without rewrites)

class SvtxEvalStack {

public:

  SvtxEvalStack(PHCompositeNode *topNode);
  virtual ~SvtxEvalStack() {}

  void next_event(PHCompositeNode *topNode);

  SvtxVertexEval*  get_vertex_eval() {return &_clustereval;}
  SvtxTrackEval*   get_track_eval() {return &_clustereval;}
  SvtxClusterEval* get_cluster_eval() {return &_clustereval;}
  SvtxHitEval*     get_hit_eval() {return _clustereval->get_rawtower_eval();}
  SvtxTruthEval*   get_truth_eval() {return _clustereval->get_truth_eval();}
  
private:
  SvtxVertexEval _vertexeval; // right now this is the top-level eval, other evals nest underneath
  SvtxTruthEval  _trutheval;  // except this one
};

#endif // __SVTXEVALSTACK_H__
