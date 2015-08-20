
#ifndef __CALORAWCLUSTEREVAL_H__
#define __CALORAWCLUSTEREVAL_H__

#include "CaloTruthEval.h"
#include "CaloRawTowerEval.h"

#include <phool/PHCompositeNode.h>
#include <g4cemc/RawCluster.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>

#include <string>
#include <set>
#include <map>

class CaloRawClusterEval {

public:

  CaloRawClusterEval(PHCompositeNode *topNode, std::string caloname);
  virtual ~CaloRawClusterEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}
  CaloTruthEval* get_truth_eval() {return _towereval.get_truth_eval();}
  CaloRawTowerEval* get_rawtower_eval() {return &_towereval;}

  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits (RawCluster* cluster);
  
  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_primaries         (RawCluster* cluster);
  PHG4Particle*           max_truth_primary_by_energy (RawCluster* cluster);

  // forwardtrace through to RawClusters
  std::set<RawCluster*> all_clusters_from(PHG4Particle* primary);
  RawCluster*           best_cluster_from(PHG4Particle* primary);
  
  // overlap calculations
  float get_energy_contribution (RawCluster* cluster, PHG4Particle* primary);
  
private:
  PHCompositeNode* _topNode;
  std::string _caloname;  
  CaloRawTowerEval _towereval;
  
  bool                                                 _do_cache;
  std::map<RawCluster*,std::set<PHG4Hit*> >            _cache_all_truth_hits;
  std::map<RawCluster*,std::set<PHG4Particle*> >       _cache_all_truth_primaries;
  std::map<RawCluster*,PHG4Particle* >                 _cache_max_truth_primary_by_energy;
  std::map<PHG4Particle*,std::set<RawCluster*> >       _cache_all_clusters_from_primary;
  std::map<PHG4Particle*,RawCluster*>                  _cache_best_cluster_from_primary;
  std::map<std::pair<RawCluster*,PHG4Particle*>,float> _cache_get_energy_contribution_primary;
};

#endif // __SVTXHITEVAL_H__
