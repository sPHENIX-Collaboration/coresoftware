#ifndef G4EVAL_SVTXCLUSTEREVAL_H
#define G4EVAL_SVTXCLUSTEREVAL_H

#include "SvtxHitEval.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>

#include <map>
#include <memory>                // for shared_ptr, less
#include <set>
#include <utility>


class PHCompositeNode;

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;

class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class SvtxTruthEval;

typedef std::multimap<float, TrkrDefs::cluskey> innerMap;

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

  // get all truth clusters matching a given layer
  std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> all_truth_clusters(TrkrDefs::cluskey cluster_key);
  
  std::pair<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> max_truth_cluster_by_energy(TrkrDefs::cluskey cluster_key);

  PHG4Hit* all_truth_hits_by_nhit(TrkrDefs::cluskey cluster);
  std::pair<int, int> gtrackid_and_layer_by_nhit(TrkrDefs::cluskey cluster);

  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles(TrkrDefs::cluskey);
  PHG4Particle* max_truth_particle_by_energy(TrkrDefs::cluskey);
  PHG4Particle* max_truth_particle_by_cluster_energy(TrkrDefs::cluskey);

  // forwardtrace through to SvtxClusters
  std::set<TrkrDefs::cluskey> all_clusters_from(PHG4Particle* truthparticle);
  std::set<TrkrDefs::cluskey> all_clusters_from(PHG4Hit* truthhit);
  TrkrDefs::cluskey best_cluster_from(PHG4Hit* truthhit);
  TrkrDefs::cluskey best_cluster_by_nhit(int gid, int layer);
  void FillRecoClusterFromG4HitCache();
  // overlap calculations
  float get_energy_contribution(TrkrDefs::cluskey cluster_key, PHG4Particle* truthparticle);
  float get_energy_contribution(TrkrDefs::cluskey cluster_key, PHG4Hit* truthhit);

  std::pair<TrkrDefs::cluskey, TrkrCluster*> reco_cluster_from_truth_cluster( TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> gclus);

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
  TrkrClusterContainer* _clustermap = nullptr;
  TrkrClusterHitAssoc* _cluster_hit_map = nullptr;
  TrkrHitTruthAssoc* _hit_truth_map = nullptr;
  PHG4TruthInfoContainer* _truthinfo = nullptr;
  PHG4HitContainer * _g4hits_tpc = nullptr;
  PHG4HitContainer * _g4hits_intt = nullptr;
  PHG4HitContainer * _g4hits_mvtx = nullptr;
  PHG4HitContainer * _g4hits_mms = nullptr;
  ActsGeometry *_tgeometry = nullptr;

  bool _strict = false ;
  int _verbosity = 0;
  unsigned int _errors = 0;

  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey cluster_key, TrkrCluster *cluster);

  bool _do_cache = true;
  std::map<TrkrDefs::cluskey, std::set<PHG4Hit*> > _cache_all_truth_hits;
  std::map<TrkrDefs::cluskey, std::map<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster> > > _cache_all_truth_clusters;
  std::map<TrkrDefs::cluskey, PHG4Hit*> _cache_max_truth_hit_by_energy;
  std::map<TrkrDefs::cluskey, std::pair<TrkrDefs::cluskey, std::shared_ptr<TrkrCluster>> > _cache_max_truth_cluster_by_energy;
  std::map<TrkrDefs::cluskey, std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<TrkrDefs::cluskey, PHG4Particle*> _cache_max_truth_particle_by_energy;
  std::map<TrkrDefs::cluskey, PHG4Particle*> _cache_max_truth_particle_by_cluster_energy;
  std::map<PHG4Particle*, std::set<TrkrDefs::cluskey> > _cache_all_clusters_from_particle;
  std::map<PHG4Hit*, std::set<TrkrDefs::cluskey> > _cache_all_clusters_from_g4hit;
  std::map<PHG4Hit*, TrkrDefs::cluskey> _cache_best_cluster_from_g4hit;
  std::map<std::pair<int, int>, TrkrDefs::cluskey> _cache_best_cluster_from_gtrackid_layer;
  std::map<std::pair<TrkrDefs::cluskey, PHG4Particle*>, float> _cache_get_energy_contribution_g4particle;
  std::map<std::pair<TrkrDefs::cluskey, PHG4Hit*>, float> _cache_get_energy_contribution_g4hit;
  std::map<std::shared_ptr<TrkrCluster>, std::pair<TrkrDefs::cluskey,TrkrCluster*> > _cache_reco_cluster_from_truth_cluster;

  // measured for low occupancy events, all in cm
  const float sig_tpc_rphi_inner = 220e-04;
  const float sig_tpc_rphi_mid = 155e-04;
  const float sig_tpc_rphi_outer = 165e-04;
  const float sig_tpc_z = 750e-04;
  const float sig_intt_rphi = 17e-04;
  const float range_intt_z = 0.9;
  const float sig_mvtx_rphi = 4.0e-04;
  const float sig_mvtx_z = 4.7e-04;
  const float sig_mms_rphi_55 = 100e-04;
  const float sig_mms_z_56 = 200e-04;


  //! cluster azimuthal searching window in _clusters_per_layer. Unit: rad
  static constexpr float _clusters_searching_window = 0.1f;
  std::multimap<unsigned int, innerMap> _clusters_per_layer;
//  std::multimap<unsigned int, PHG4Hit*> _g4hits_per_layer;
};

#endif  // G4EVAL_SVTXCLUSTEREVAL_H
