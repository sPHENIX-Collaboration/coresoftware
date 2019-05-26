#ifndef G4EVAL_CALORAWCLUSTEREVAL_H
#define G4EVAL_CALORAWCLUSTEREVAL_H

#include "CaloRawTowerEval.h"
#include "CaloTruthEval.h"

#include <map>
#include <set>
#include <string>
#include <utility>

class PHCompositeNode;
class PHG4Hit;
class PHG4Particle;
class PHG4Shower;
class RawClusterContainer;
class RawCluster;
class RawTowerContainer;

class CaloRawClusterEval
{
 public:
  /// example caloname: CEMC, HCALIN, HCALOUT
  CaloRawClusterEval(PHCompositeNode* topNode, const std::string& caloname);
  virtual ~CaloRawClusterEval();

  /// get the hash id for this calorimeter volume
  int get_caloid() { return get_truth_eval()->get_caloid(); }

  /// reinitialize the eval for a new event
  void next_event(PHCompositeNode* topNode);

  /// activate or deactivate the memory caching inside the evaluation module
  void do_caching(bool do_cache)
  {
    _do_cache = do_cache;
    _towereval.do_caching(do_cache);
  }

  /// strict mode will assert when an error is detected
  /// non-strict mode will notice and report at the End()
  void set_strict(bool strict)
  {
    _strict = strict;
    _towereval.set_strict(strict);
  }

  /// get a count of the errors discovered thus far
  unsigned int get_errors() { return _errors + _towereval.get_errors(); }

  /// adjust the messaging from the evalutaion module
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _towereval.set_verbosity(verbosity);
  }

  /// get a copy of the lower level truth eval and its memory cache
  CaloTruthEval* get_truth_eval() { return _towereval.get_truth_eval(); }

  /// get a copy of the lower level tower eval and its memory cache
  CaloRawTowerEval* get_rawtower_eval() { return &_towereval; }

  // ---reduced sim node or better----------------------------------------------

  /// has the eval initialized correctly for reduced sim DST nodes?
  bool has_reduced_node_pointers();

  // shower interface

  /// what primary showers contributed energy to this cluster?
  std::set<PHG4Shower*> all_truth_primary_showers(RawCluster* cluster);

  /// which primary shower contributed the most energy to this cluster?
  PHG4Shower* max_truth_primary_shower_by_energy(RawCluster* cluster);

  /// what clusters did this primary truth shower contribute energy to?
  std::set<RawCluster*> all_clusters_from(PHG4Shower* primary);

  /// which cluster did this primary truth shower contribute the most energy to?
  RawCluster* best_cluster_from(PHG4Shower* primary);

  /// how much energy did this primary truth shower contribute to this cluster
  float get_energy_contribution(RawCluster* cluster, PHG4Shower* primary);

  // particle interface

  /// what particles contributed energy to this cluster?
  std::set<PHG4Particle*> all_truth_primary_particles(RawCluster* cluster);

  /// which particle contributed the most energy to this cluster?
  PHG4Particle* max_truth_primary_particle_by_energy(RawCluster* cluster);

  /// what clusters did this primary truth particle contribute energy to?
  std::set<RawCluster*> all_clusters_from(PHG4Particle* primary);

  /// which cluster did this primary truth particle contribute the most energy to?
  RawCluster* best_cluster_from(PHG4Particle* primary);

  /// how much energy did this primary truth particle contribute to this cluster
  float get_energy_contribution(RawCluster* cluster, PHG4Particle* primary);

  // ---full sim node required--------------------------------------------------

  /// has the eval initialized correctly for full sim DST nodes?
  bool has_full_node_pointers();

  /// what truth hits contributed energy to this tower?
  std::set<PHG4Hit*> all_truth_hits(RawCluster* cluster);

 private:
  void get_node_pointers(PHCompositeNode* topNode);

  std::string _caloname;
  CaloRawTowerEval _towereval;
  RawClusterContainer* _clusters;
  RawTowerContainer* _towers;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::map<RawCluster*, std::set<PHG4Shower*> > _cache_all_truth_primary_showers;
  std::map<RawCluster*, PHG4Shower*> _cache_max_truth_primary_shower_by_energy;
  std::map<PHG4Shower*, std::set<RawCluster*> > _cache_all_clusters_from_primary_shower;
  std::map<PHG4Shower*, RawCluster*> _cache_best_cluster_from_primary_shower;
  std::map<std::pair<RawCluster*, PHG4Shower*>, float> _cache_get_energy_contribution_primary_shower;

  std::map<RawCluster*, std::set<PHG4Particle*> > _cache_all_truth_primary_particles;
  std::map<RawCluster*, PHG4Particle*> _cache_max_truth_primary_particle_by_energy;
  std::map<PHG4Particle*, std::set<RawCluster*> > _cache_all_clusters_from_primary_particle;
  std::map<PHG4Particle*, RawCluster*> _cache_best_cluster_from_primary_particle;
  std::map<std::pair<RawCluster*, PHG4Particle*>, float> _cache_get_energy_contribution_primary_particle;

  std::map<RawCluster*, std::set<PHG4Hit*> > _cache_all_truth_hits;
};

#endif  // G4EVAL_CALORAWCLUSTEREVAL_H
