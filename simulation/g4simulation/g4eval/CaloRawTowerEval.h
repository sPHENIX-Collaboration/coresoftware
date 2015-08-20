
#ifndef __CALORAWTOWEREVAL_H__
#define __CALORAWTOWEREVAL_H__

#include "CaloTruthEval.h"

#include <phool/PHCompositeNode.h>
#include <g4cemc/RawTower.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>

#include <string>
#include <set>
#include <map>

class CaloRawTowerEval {

public:

  CaloRawTowerEval(PHCompositeNode *topNode, std::string caloname);
  virtual ~CaloRawTowerEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}
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
  float get_energy_contribution (RawTower* svtxhit, PHG4Particle* primary);
  
private:
  PHCompositeNode* _topNode;
  std::string _caloname;
  CaloTruthEval _trutheval;

  bool                                               _do_cache;
  std::map<RawTower*,std::set<PHG4Hit*> >            _cache_all_truth_hits;
  std::map<RawTower*,std::set<PHG4Particle*> >       _cache_all_truth_primaries;
  std::map<RawTower*,PHG4Particle* >                 _cache_max_truth_primary_by_energy;
  std::map<PHG4Particle*,std::set<RawTower*> >       _cache_all_towers_from_primary;
  std::map<PHG4Particle*,RawTower*>                  _cache_best_tower_from_primary;
  std::map<std::pair<RawTower*,PHG4Particle*>,float> _cache_get_energy_contribution_primary;
};

#endif // __SVTXHITEVAL_H__
