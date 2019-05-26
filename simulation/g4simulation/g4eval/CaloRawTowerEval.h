#ifndef G4EVAL_CALORAWTOWEREVAL_H
#define G4EVAL_CALORAWTOWEREVAL_H

#include "CaloTruthEval.h"

#include <map>
#include <set>
#include <string>
#include <utility>

class PHCompositeNode;

class PHG4CellContainer;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4Shower;
class PHG4TruthInfoContainer;

class RawTower;
class RawTowerContainer;

class CaloRawTowerEval
{
 public:
  /// example caloname: CEMC, HCALIN, HCALOUT
  CaloRawTowerEval(PHCompositeNode* topNode, const std::string& caloname);
  virtual ~CaloRawTowerEval();

  /// get the hash id for this calorimeter volume
  int get_caloid() { return get_truth_eval()->get_caloid(); }

  /// reinitialize the eval for a new event
  void next_event(PHCompositeNode* topNode);

  /// activate or deactivate the memory caching inside the evaluation module
  void do_caching(bool do_cache)
  {
    _do_cache = do_cache;
    _trutheval.do_caching(do_cache);
  }

  /// strict mode will assert when an error is detected
  /// non-strict mode will notice and report at the End()
  void set_strict(bool strict)
  {
    _strict = strict;
    _trutheval.set_strict(strict);
  }

  /// get a count of the errors discovered thus far
  unsigned int get_errors() { return _errors + _trutheval.get_errors(); }

  /// adjust the messaging from the evalutaion module
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _trutheval.set_verbosity(verbosity);
  }

  /// get a copy of the lower level truth eval and its memory cache
  CaloTruthEval* get_truth_eval() { return &_trutheval; }

  // ---reduced sim node or better----------------------------------------------

  /// has the eval initialized correctly for reduced sim DST nodes?
  bool has_reduced_node_pointers();

  // shower interface

  /// what primary showers contributed energy to this tower?
  std::set<PHG4Shower*> all_truth_primary_showers(RawTower* tower);

  /// which primary shower contributed the most energy to this tower?
  PHG4Shower* max_truth_primary_shower_by_energy(RawTower* tower);

  /// what towers did this primary shower contribute energy to?
  std::set<RawTower*> all_towers_from(PHG4Shower* primary);

  /// which tower did this primary shower contribute the most energy to?
  RawTower* best_tower_from(PHG4Shower* primary);

  /// how much energy did this primary shower contribute to this tower?
  float get_energy_contribution(RawTower* tower, PHG4Shower* primary);

  // particle interface

  /// what particles contributed energy to this tower?
  std::set<PHG4Particle*> all_truth_primary_particles(RawTower* tower);

  /// which particle contributed the most energy to this tower?
  PHG4Particle* max_truth_primary_particle_by_energy(RawTower* tower);

  /// what towers did this primary truth particle contribute energy to?
  std::set<RawTower*> all_towers_from(PHG4Particle* primary);

  /// which tower did the primary truth particle contribute the most energy to?
  RawTower* best_tower_from(PHG4Particle* primary);

  /// how much energy did this primary truth particle contribute to this tower?
  float get_energy_contribution(RawTower* tower, PHG4Particle* primary);

  // ---full sim node required--------------------------------------------------

  /// has the eval initialized correctly for full sim DST nodes?
  bool has_full_node_pointers();

  /// what truth hits contributed energy to this tower?
  std::set<PHG4Hit*> all_truth_hits(RawTower* tower);

 private:
  void get_node_pointers(PHCompositeNode* topNode);

  std::string _caloname;
  CaloTruthEval _trutheval;
  RawTowerContainer* _towers;
  PHG4CellContainer* _g4cells;
  PHG4HitContainer* _g4hits;
  PHG4TruthInfoContainer* _truthinfo;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;

  std::map<RawTower*, std::set<PHG4Shower*> > _cache_all_truth_primary_showers;
  std::map<RawTower*, PHG4Shower*> _cache_max_truth_primary_shower_by_energy;
  std::map<PHG4Shower*, std::set<RawTower*> > _cache_all_towers_from_primary_shower;
  std::map<PHG4Shower*, RawTower*> _cache_best_tower_from_primary_shower;
  std::map<std::pair<RawTower*, PHG4Shower*>, float> _cache_get_energy_contribution_primary_shower;

  std::map<RawTower*, std::set<PHG4Particle*> > _cache_all_truth_primary_particles;
  std::map<RawTower*, PHG4Particle*> _cache_max_truth_primary_particle_by_energy;
  std::map<PHG4Particle*, std::set<RawTower*> > _cache_all_towers_from_primary_particle;
  std::map<PHG4Particle*, RawTower*> _cache_best_tower_from_primary_particle;
  std::map<std::pair<RawTower*, PHG4Particle*>, float> _cache_get_energy_contribution_primary_particle;

  std::map<RawTower*, std::set<PHG4Hit*> > _cache_all_truth_hits;
};

#endif  // G4EVAL_CALORAWTOWEREVAL_H
