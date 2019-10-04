#ifndef G4EVAL_SVTXCLUSTEREVAL_H
#define G4EVAL_SVTXCLUSTEREVAL_H

#include "SvtxHitEval.h"

#include <trackbase/TrkrDefs.h>

#include <map>
#include <set>
#include <utility>


class PHCompositeNode;

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;

class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class SvtxTruthEval;

using namespace std;
typedef multimap<float, TrkrDefs::cluskey> innerMap;

class SvtxClusterEval
{
 public:
  SvtxClusterEval(PHCompositeNode* topNode);
  virtual ~SvtxClusterEval();

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache)
  {
    _do_cache = do_cache;
    _hiteval.do_caching(do_cache);
  }
  void set_strict(bool strict)
  {
    _strict = strict;
    _hiteval.set_strict(strict);
  }
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _hiteval.set_verbosity(verbosity);
  }

  // access the clustereval (and its cached values)
  SvtxHitEval* get_hit_eval() { return &_hiteval; }
  SvtxTruthEval* get_truth_eval() { return _hiteval.get_truth_eval(); }

  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits(TrkrDefs::cluskey cluster);
  PHG4Hit* max_truth_hit_by_energy(TrkrDefs::cluskey);

  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles(TrkrDefs::cluskey);
  PHG4Particle* max_truth_particle_by_energy(TrkrDefs::cluskey);

  // forwardtrace through to SvtxClusters
  std::set<TrkrDefs::cluskey> all_clusters_from(PHG4Particle* truthparticle);
  std::set<TrkrDefs::cluskey> all_clusters_from(PHG4Hit* truthhit);
  TrkrDefs::cluskey best_cluster_from(PHG4Hit* truthhit);

  // overlap calculations
  float get_energy_contribution(TrkrDefs::cluskey cluster_key, PHG4Particle* truthparticle);
  float get_energy_contribution(TrkrDefs::cluskey cluster_key, PHG4Hit* truthhit);

  unsigned int get_errors() { return _errors + _hiteval.get_errors(); }

 private:
  void get_node_pointers(PHCompositeNode* topNode);
  void fill_cluster_layer_map();
  //  void fill_g4hit_layer_map();
  bool has_node_pointers();

  //! Fast approximation of atan2() for cluster searching
  //! From https://www.dsprelated.com/showarticle/1052.php
  float fast_approx_atan2(float y, float x);
  float fast_approx_atan2(float y2x);

  SvtxHitEval _hiteval;
  TrkrClusterContainer* _clustermap;
  TrkrClusterHitAssoc* _cluster_hit_map;
  TrkrHitTruthAssoc* _hit_truth_map;
  PHG4TruthInfoContainer* _truthinfo;
  PHG4HitContainer * _g4hits_tpc;
  PHG4HitContainer * _g4hits_intt;
  PHG4HitContainer * _g4hits_mvtx;
  PHG4HitContainer * _g4hits_outertracker;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::map<TrkrDefs::cluskey, std::set<PHG4Hit*> > _cache_all_truth_hits;
  std::map<TrkrDefs::cluskey, PHG4Hit*> _cache_max_truth_hit_by_energy;
  std::map<TrkrDefs::cluskey, std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<TrkrDefs::cluskey, PHG4Particle*> _cache_max_truth_particle_by_energy;
  std::map<PHG4Particle*, std::set<TrkrDefs::cluskey> > _cache_all_clusters_from_particle;
  std::map<PHG4Hit*, std::set<TrkrDefs::cluskey> > _cache_all_clusters_from_g4hit;
  std::map<PHG4Hit*, TrkrDefs::cluskey> _cache_best_cluster_from_g4hit;
  std::map<std::pair<TrkrDefs::cluskey, PHG4Particle*>, float> _cache_get_energy_contribution_g4particle;
  std::map<std::pair<TrkrDefs::cluskey, PHG4Hit*>, float> _cache_get_energy_contribution_g4hit;

#if !defined(__CINT__) || defined(__CLING__)
  //! cluster azimuthal searching window in _clusters_per_layer. Unit: rad
  static constexpr float _clusters_searching_window = 0.1f;
  std::multimap<unsigned int, innerMap> _clusters_per_layer;
//  std::multimap<unsigned int, PHG4Hit*> _g4hits_per_layer;
#endif
};

#endif  // G4EVAL_SVTXCLUSTEREVAL_H
