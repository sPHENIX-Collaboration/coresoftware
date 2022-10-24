#ifndef G4EVAL_SVTXHITEVAL_H
#define G4EVAL_SVTXHITEVAL_H

#include "SvtxTruthEval.h"

#include <trackbase/TrkrDefs.h>

#include <map>
#include <set>
#include <utility>

class PHCompositeNode;

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;

class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrHitTruthAssoc;

class SvtxHitEval
{
 public:
  SvtxHitEval(PHCompositeNode* topNode);
  virtual ~SvtxHitEval();

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache)
  {
    _do_cache = do_cache;
    _trutheval.do_caching(do_cache);
  }
  void set_strict(bool strict)
  {
    _strict = strict;
    _trutheval.set_strict(strict);
  }
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _trutheval.set_verbosity(verbosity);
  }

  // access the clustereval (and its cached values)
  SvtxTruthEval* get_truth_eval() { return &_trutheval; }

  //PHG4Cell* get_cell(SvtxHit* hit);

  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits(TrkrDefs::hitkey hit_key);
  PHG4Hit* max_truth_hit_by_energy(TrkrDefs::hitkey hit_key);

  // backtrace through to PHG4Hits for a specific tracker
  std::set<PHG4Hit*> all_truth_hits(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid);
  PHG4Hit* max_truth_hit_by_energy(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid);

  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles(TrkrDefs::hitkey hit_key);
  PHG4Particle* max_truth_particle_by_energy(TrkrDefs::hitkey hit_key);

  // backtrace through to PHG4Particles for a specific tracker
  std::set<PHG4Particle*> all_truth_particles(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid);
  PHG4Particle* max_truth_particle_by_energy(TrkrDefs::hitkey hit_key, const TrkrDefs::TrkrId trkrid);

  // forwardtrace through to SvtxHits
  std::set<TrkrDefs::hitkey> all_hits_from(PHG4Particle* truthparticle);
  std::set<TrkrDefs::hitkey> all_hits_from(PHG4Hit* truthhit);
  TrkrDefs::hitkey best_hit_from(PHG4Hit* truthhit);

  // overlap calculations
  float get_energy_contribution(TrkrDefs::hitkey, PHG4Particle* truthparticle);
  float get_energy_contribution(TrkrDefs::hitkey, PHG4Hit* truthhit);

  unsigned int get_errors() { return _errors + _trutheval.get_errors(); }

 private:
  void get_node_pointers(PHCompositeNode* topNode);
  bool has_node_pointers();

  SvtxTruthEval _trutheval;
  TrkrHitSetContainer* _hitmap;
  TrkrClusterContainer* _clustermap;
  TrkrHitTruthAssoc* _hit_truth_map;

  PHG4HitContainer* _g4hits_tpc;
  PHG4HitContainer* _g4hits_intt;
  PHG4HitContainer* _g4hits_mvtx;
  PHG4HitContainer* _g4hits_mms;

  PHG4TruthInfoContainer* _truthinfo;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::map<TrkrDefs::hitkey, std::set<PHG4Hit*> > _cache_all_truth_hits;
  std::map<TrkrDefs::hitkey, PHG4Hit*> _cache_max_truth_hit_by_energy;
  std::map<TrkrDefs::hitkey, std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<TrkrDefs::hitkey, PHG4Particle*> _cache_max_truth_particle_by_energy;
  std::map<PHG4Particle*, std::set<TrkrDefs::hitkey> > _cache_all_hits_from_particle;
  std::map<PHG4Hit*, std::set<TrkrDefs::hitkey> > _cache_all_hits_from_g4hit;
  std::map<PHG4Hit*, TrkrDefs::hitkey> _cache_best_hit_from_g4hit;
  std::map<std::pair<TrkrDefs::hitkey, PHG4Particle*>, float> _cache_get_energy_contribution_g4particle;
  std::map<std::pair<TrkrDefs::hitkey, PHG4Hit*>, float> _cache_get_energy_contribution_g4hit;
};

#endif  // G4EVAL_SVTXHITEVAL_H
