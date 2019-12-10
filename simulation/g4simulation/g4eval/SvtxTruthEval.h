#ifndef G4EVAL_SVTXTRUTHEVAL_H
#define G4EVAL_SVTXTRUTHEVAL_H

#include "BaseTruthEval.h"

class PHCompositeNode;

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;

#include <map>
#include <set>

class SvtxTruthEval
{
 public:
  SvtxTruthEval(PHCompositeNode* topNode);
  virtual ~SvtxTruthEval();

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache) { _do_cache = do_cache; }
  void set_strict(bool strict)
  {
    _strict = strict;
    _basetrutheval.set_strict(strict);
  }
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _basetrutheval.set_verbosity(verbosity);
  }

  std::set<PHG4Hit*> all_truth_hits();
  std::set<PHG4Hit*> all_truth_hits(PHG4Particle* particle);
  PHG4Particle* get_particle(PHG4Hit* g4hit);
  int get_embed(PHG4Particle* particle);
  PHG4VtxPoint* get_vertex(PHG4Particle* particle);
  bool is_primary(PHG4Particle* particle);
  PHG4Particle* get_primary_particle(PHG4Hit* g4hit);
  PHG4Particle* get_primary_particle(PHG4Particle* particle);

  bool is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle);
  bool are_same_particle(PHG4Particle* p1, PHG4Particle* p2);
  bool are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2);

  PHG4Hit* get_innermost_truth_hit(PHG4Particle* particle);
  PHG4Hit* get_outermost_truth_hit(PHG4Particle* particle);

  unsigned int get_errors() { return _errors + _basetrutheval.get_errors(); }

 private:
  void get_node_pointers(PHCompositeNode* topNode);
  bool has_node_pointers();

  BaseTruthEval _basetrutheval;

  PHG4TruthInfoContainer* _truthinfo;
  PHG4HitContainer* _g4hits_svtx;
  PHG4HitContainer* _g4hits_tracker;
  PHG4HitContainer* _g4hits_maps;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::set<PHG4Hit*> _cache_all_truth_hits;
  std::map<PHG4Particle*, std::set<PHG4Hit*> > _cache_all_truth_hits_g4particle;
  std::map<PHG4Particle*, PHG4Hit*> _cache_get_innermost_truth_hit;
  std::map<PHG4Particle*, PHG4Hit*> _cache_get_outermost_truth_hit;
  std::map<PHG4Hit*, PHG4Particle*> _cache_get_primary_particle_g4hit;
};

#endif  // G4EVAL_SVTXTRUTHEVAL_H
