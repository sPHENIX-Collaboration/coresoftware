#ifndef G4EVAL_SVTXVERTEXEVAL_H
#define G4EVAL_SVTXVERTEXEVAL_H

#include "SvtxTrackEval.h"

#include <map>
#include <set>
#include <string>  // for string
#include <utility>

class PHCompositeNode;

class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;

class SvtxClusterEval;
class SvtxHitEval;
class SvtxTrackMap;
class SvtxTruthEval;
class SvtxVertex;
class SvtxVertexMap;

class SvtxVertexEval
{
 public:
  SvtxVertexEval(PHCompositeNode* topNode);
  virtual ~SvtxVertexEval();

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache)
  {
    _do_cache = do_cache;
    _trackeval.do_caching(do_cache);
  }
  void set_strict(bool strict)
  {
    _strict = strict;
    _trackeval.set_strict(strict);
  }
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _trackeval.set_verbosity(verbosity);
  }

  // access the sub evals (and the cached values)
  SvtxTrackEval* get_track_eval() { return &_trackeval; }
  SvtxClusterEval* get_cluster_eval() { return _trackeval.get_cluster_eval(); }
  SvtxHitEval* get_hit_eval() { return _trackeval.get_hit_eval(); }
  SvtxTruthEval* get_truth_eval() { return _trackeval.get_truth_eval(); }

  // backtrace through to PHG4Hits
  std::set<PHG4Particle*> all_truth_particles(SvtxVertex* vertex);

  // backtrace through to PHG4VtxPoints
  std::set<PHG4VtxPoint*> all_truth_points(SvtxVertex* vertex);
  PHG4VtxPoint* max_truth_point_by_ntracks(SvtxVertex* vertex);

  // forwardtrace through to SvtxVertexs
  std::set<SvtxVertex*> all_vertexes_from(PHG4VtxPoint* truthpoint);
  SvtxVertex* best_vertex_from(PHG4VtxPoint* truthpoint);

  // overlap calculations
  unsigned int get_ntracks_contribution(SvtxVertex* svtxvertex, PHG4VtxPoint* truthpoint);

  unsigned int get_errors() { return _errors + _trackeval.get_errors(); }

  void set_use_initial_vertex(bool use_init_vertex) { _use_initial_vertex = use_init_vertex; }
  void set_use_genfit_vertex(bool use_genfit_vertex) { _use_genfit_vertex = use_genfit_vertex; }

  void set_track_nodename(const std::string& name);

 private:
  void get_node_pointers(PHCompositeNode* topNode);
  bool has_node_pointers();

  SvtxTrackEval _trackeval;
  SvtxVertexMap* _vertexmap;
  SvtxTrackMap* _trackmap;
  PHG4TruthInfoContainer* _truthinfo;

  bool _strict;
  bool _use_initial_vertex = true;
  bool _use_genfit_vertex = false;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::map<SvtxVertex*, std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<SvtxVertex*, std::set<PHG4VtxPoint*> > _cache_all_truth_points;
  std::map<SvtxVertex*, PHG4VtxPoint*> _cache_max_truth_point_by_ntracks;
  std::map<PHG4VtxPoint*, std::set<SvtxVertex*> > _cache_all_vertexes_from_point;
  std::map<PHG4VtxPoint*, SvtxVertex*> _cache_best_vertex_from_point;
  std::map<std::pair<SvtxVertex*, PHG4VtxPoint*>, unsigned int> _cache_get_ntracks_contribution;
  std::string m_TrackNodeName;
};

#endif  // G4EVAL_SVTXVERTEXEVAL_H
