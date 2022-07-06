#ifndef G4EVAL_SVTXEVALSTACK_H
#define G4EVAL_SVTXEVALSTACK_H

#include "SvtxVertexEval.h"

class PHCompositeNode;
class SvtxClusterEval;
class SvtxHitEval;
class SvtxTrackEval;
class SvtxTruthEval;

#include <string>  // for string

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
  //void set_over_write_vertexmap(bool over_write) {_vertexeval.set_over_write_vertexmap(over_write);}
  void set_use_initial_vertex(bool use_init_vtx) { _vertexeval.set_use_initial_vertex(use_init_vtx); }
  void set_use_genfit_vertex(bool use_genfit_vtx) { _vertexeval.set_use_genfit_vertex(use_genfit_vtx); }
  void set_verbosity(int verbosity) { _vertexeval.set_verbosity(verbosity); }

  SvtxVertexEval* get_vertex_eval() { return &_vertexeval; }
  SvtxTrackEval* get_track_eval() { return _vertexeval.get_track_eval(); }
  SvtxClusterEval* get_cluster_eval() { return _vertexeval.get_cluster_eval(); }
  SvtxHitEval* get_hit_eval() { return _vertexeval.get_hit_eval(); }
  SvtxTruthEval* get_truth_eval() { return _vertexeval.get_truth_eval(); }

  unsigned int get_errors() { return _vertexeval.get_errors(); }

  void set_track_nodename(const std::string& name) { _vertexeval.set_track_nodename(name); }

 private:
  SvtxVertexEval _vertexeval;  // right now this is the top-level eval

};

#endif  // G4EVAL_SVTXEVALSTACK_H
