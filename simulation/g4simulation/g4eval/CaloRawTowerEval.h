
#ifndef __CALORAWTOWEREVAL_H__
#define __CALORAWTOWEREVAL_H__

#include "CaloTruthEval.h"

#include <phool/PHCompositeNode.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <string>
#include <set>
#include <map>

class CaloRawTowerEval {

public:

  CaloRawTowerEval(PHCompositeNode *topNode, std::string caloname);
  virtual ~CaloRawTowerEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {
    _do_cache = do_cache;
    _trutheval.do_caching(do_cache);
  }
  void set_strict(bool strict) {
    _strict = strict;
    _trutheval.set_strict(strict);
  }

  CaloTruthEval* get_truth_eval() {return &_trutheval;}

  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits (RawTower* tower);
  
  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_primaries         (RawTower* tower);
  PHG4Particle*           max_truth_primary_by_energy (RawTower* tower);

  // forwardtrace through to RawTowers
  std::set<RawTower*> all_towers_from(PHG4Particle* primary);
  RawTower*           best_tower_from(PHG4Particle* primary);
  
  // overlap calculations
  float get_energy_contribution (RawTower* tower, PHG4Particle* primary);
  
private:

  void get_node_pointers(PHCompositeNode *topNode);

  std::string _caloname;
  CaloTruthEval _trutheval;  
  RawTowerContainer* _towers;
  PHG4CylinderCellContainer* _g4cells;
  PHG4HitContainer* _g4hits;
  PHG4TruthInfoContainer* _truthinfo;

  bool _strict;
  
  bool                                               _do_cache;
  std::map<RawTower*,std::set<PHG4Hit*> >            _cache_all_truth_hits;
  std::map<RawTower*,std::set<PHG4Particle*> >       _cache_all_truth_primaries;
  std::map<RawTower*,PHG4Particle* >                 _cache_max_truth_primary_by_energy;
  std::map<PHG4Particle*,std::set<RawTower*> >       _cache_all_towers_from_primary;
  std::map<PHG4Particle*,RawTower*>                  _cache_best_tower_from_primary;
  std::map<std::pair<RawTower*,PHG4Particle*>,float> _cache_get_energy_contribution_primary;
};

#endif // __SVTXHITEVAL_H__
