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

class Vertex;
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
  std::set<PHG4Particle*> all_truth_particles(const Vertex* vertex);

  // backtrace through to PHG4VtxPoints
  std::set<PHG4VtxPoint*> all_truth_points(const Vertex* vertex);
  PHG4VtxPoint* max_truth_point_by_ntracks(const Vertex* vertex);

  // forwardtrace through to SvtxVertexs
  std::set<const Vertex*> all_vertexes_from(PHG4VtxPoint* truthpoint);
  const Vertex* best_vertex_from(PHG4VtxPoint* truthpoint);

  // overlap calculations
  unsigned int get_ntracks_contribution(const Vertex* svtxvertex, PHG4VtxPoint* truthpoint);

  unsigned int get_errors() { return _errors + _trackeval.get_errors(); }

  void set_use_initial_vertex(bool use_init_vertex) { _use_initial_vertex = use_init_vertex; }
  void set_use_genfit_vertex(bool use_genfit_vertex) { _use_genfit_vertex = use_genfit_vertex; }

  void set_track_nodename(const std::string& name);

 private:
  void get_node_pointers(PHCompositeNode* topNode);
  bool has_node_pointers();

  SvtxTrackEval _trackeval;
  SvtxVertexMap* _vertexmap = nullptr;
  SvtxTrackMap* _trackmap = nullptr;
  PHG4TruthInfoContainer* _truthinfo = nullptr;

  bool _strict = false;
  bool _use_initial_vertex = true;
  bool _use_genfit_vertex = false;
  int _verbosity = 0;
  unsigned int _errors = 0;

  bool _do_cache = true;
  std::map<const Vertex*, std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<const Vertex*, std::set<PHG4VtxPoint*> > _cache_all_truth_points;
  std::map<const Vertex*, PHG4VtxPoint*> _cache_max_truth_point_by_ntracks;
  std::map<PHG4VtxPoint*, std::set<const Vertex*> > _cache_all_vertexes_from_point;
  std::map<PHG4VtxPoint*, const Vertex*> _cache_best_vertex_from_point;
  std::map<std::pair<const Vertex*, PHG4VtxPoint*>, unsigned int> _cache_get_ntracks_contribution;
  std::string m_TrackNodeName;
};

#endif  // G4EVAL_SVTXVERTEXEVAL_H
