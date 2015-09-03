
#ifndef __JETRECOEVAL_H__
#define __JETRECOEVAL_H__

#include "CaloTruthEval.h"
#include "CaloRawTowerEval.h"

#include <phool/PHCompositeNode.h>
#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>

#include <string>
#include <set>
#include <map>

class JetrecoEval {

public:

  JetRecoEval(PHCompositeNode *topNode, std::string jetname);
  virtual ~JetRecoEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}

  JetTruthEval*     get_truth_eval() {return &_jettrutheval;}
  SvtxEvalStack*    get_svtx_eval_stack();
  CaloEvalStack*    get_calo_eval_stack(std::string caloname) {return &_caloevalstacks[caloname];}

  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits (Jet* recojet);

  // backtrace through to truth jets
  std::set<Jet*> all_truth_jets          (Jet* recojet);
  Jet*           max_truth_jet_by_energy (Jet* recojet);
  
  // backtrace through to truth particles
  std::set<PHG4Particle*> all_truth_particles (Jet* recojet);

  // forwardtrace through to Reco Jets
  std::set<Jet*> all_jets_from(Jet* truthjet);
  Jet*           best_jet_from(Jet* truthjet);
  
  // overlap calculations (to reco from truth)
  float get_energy_contribution (Jet* recojet, Jet* truthjet);
  
private:

  void get_node_pointers(PHCompositeNode *topNode);
  
  std::string _jetname;  

  JetTruthEval _jettrutheval;
  SvtxEvalStack _svtxevalstack;
  std::map<std::string,CaloEvalStack> _caloevalstacks;
  
  bool                                                 _do_cache;
  std::map<RawCluster*,std::set<PHG4Hit*> >            _cache_all_truth_hits;
  std::map<RawCluster*,std::set<PHG4Particle*> >       _cache_all_truth_primaries;
  std::map<RawCluster*,PHG4Particle* >                 _cache_max_truth_primary_by_energy;
  std::map<PHG4Particle*,std::set<RawCluster*> >       _cache_all_clusters_from_primary;
  std::map<PHG4Particle*,RawCluster*>                  _cache_best_cluster_from_primary;
  std::map<std::pair<RawCluster*,PHG4Particle*>,float> _cache_get_energy_contribution_primary;
};

#endif // __SVTXHITEVAL_H__
