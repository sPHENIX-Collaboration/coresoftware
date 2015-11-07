
#ifndef __SVTXCLUSTEREVAL_H__
#define __SVTXCLUSTEREVAL_H__

#include "SvtxHitEval.h"
#include "SvtxTruthEval.h"

#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxHitMap.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <set>
#include <map>

class SvtxClusterEval {

public:

  SvtxClusterEval(PHCompositeNode *topNode);
  virtual ~SvtxClusterEval();

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {
    _do_cache = do_cache;
    _hiteval.do_caching(do_cache);
  }
  void set_strict(bool strict) {
    _strict = strict;
    _hiteval.set_strict(strict);
  }
  void set_verbosity(int verbosity) {
    _verbosity = verbosity;
    _hiteval.set_verbosity(verbosity);
  }
  
  // access the clustereval (and its cached values)
  SvtxHitEval* get_hit_eval() {return &_hiteval;}
  SvtxTruthEval* get_truth_eval() {return _hiteval.get_truth_eval();}
  
  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits          (SvtxCluster* cluster);
  PHG4Hit*           max_truth_hit_by_energy (SvtxCluster* cluster);
  
  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles          (SvtxCluster* cluster);
  PHG4Particle*           max_truth_particle_by_energy (SvtxCluster* cluster);

  // forwardtrace through to SvtxClusters
  std::set<SvtxCluster*> all_clusters_from(PHG4Particle* truthparticle);  
  std::set<SvtxCluster*> all_clusters_from(PHG4Hit* truthhit);
  SvtxCluster*           best_cluster_from(PHG4Hit* truthhit);
  
  // overlap calculations
  float get_energy_contribution (SvtxCluster* svtxcluster, PHG4Particle* truthparticle);
  float get_energy_contribution (SvtxCluster* svtxcluster, PHG4Hit* truthhit);

  unsigned int get_errors() {return _errors + _hiteval.get_errors();}
  
private:

  void get_node_pointers(PHCompositeNode* topNode);
  
  SvtxHitEval _hiteval;
  SvtxClusterMap* _clustermap;
  SvtxHitMap* _hitmap;
  PHG4TruthInfoContainer* _truthinfo;

  bool _strict;
  int _verbosity;
  unsigned int _errors;
  
  bool                                                  _do_cache;
  std::map<SvtxCluster*,std::set<PHG4Hit*> >            _cache_all_truth_hits;
  std::map<SvtxCluster*,PHG4Hit*>                       _cache_max_truth_hit_by_energy;
  std::map<SvtxCluster*,std::set<PHG4Particle*> >       _cache_all_truth_particles;
  std::map<SvtxCluster*,PHG4Particle* >                 _cache_max_truth_particle_by_energy;
  std::map<PHG4Particle*,std::set<SvtxCluster*> >       _cache_all_clusters_from_particle;
  std::map<PHG4Hit*,std::set<SvtxCluster*> >            _cache_all_clusters_from_g4hit;
  std::map<PHG4Hit*,SvtxCluster* >                      _cache_best_cluster_from_g4hit;
  std::map<std::pair<SvtxCluster*,PHG4Particle*>,float> _cache_get_energy_contribution_g4particle;
  std::map<std::pair<SvtxCluster*,PHG4Hit*>,float>      _cache_get_energy_contribution_g4hit;
};

#endif // __SVTXCLUSTEREVAL_H__
