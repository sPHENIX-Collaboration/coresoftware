#ifndef G4EVAL_SVTXEVALSTACK_H
#define G4EVAL_SVTXEVALSTACK_H

#include "SvtxVertexEval.h"

class PHCompositeNode;
class SvtxClusterEval;
class SvtxHitEval;
class SvtxTrackEval;
class SvtxTruthEval;

// This user class provides pointers to the
// full set of tracking evaluators and
// protects the user from future introduction
// of new eval heirachies (new eval objects can
// be introduced without rewrites)

class SvtxEvalStack
{
 public:
  SvtxEvalStack(PHCompositeNode* topNode);
  virtual ~SvtxEvalStack() {}

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache) { _vertexeval.do_caching(do_cache); }
  void set_strict(bool strict) { _vertexeval.set_strict(strict); }
  void set_verbosity(int verbosity) { _vertexeval.set_verbosity(verbosity); }

  SvtxVertexEval* get_vertex_eval() { return &_vertexeval; }
  SvtxTrackEval* get_track_eval() { return _vertexeval.get_track_eval(); }
  SvtxClusterEval* get_cluster_eval() { return _vertexeval.get_cluster_eval(); }
  SvtxHitEval* get_hit_eval() { return _vertexeval.get_hit_eval(); }
  SvtxTruthEval* get_truth_eval() { return _vertexeval.get_truth_eval(); }

  unsigned int get_errors() { return _vertexeval.get_errors(); }

 private:
  SvtxVertexEval _vertexeval;  // right now this is the top-level eval
};

#endif  // G4EVAL_SVTXEVALSTACK_H
