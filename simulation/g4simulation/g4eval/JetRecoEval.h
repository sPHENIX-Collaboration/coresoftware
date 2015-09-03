
#ifndef __JETRECOEVAL_H__
#define __JETRECOEVAL_H__

#include "JetTruthEval.h"

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

class JetRecoEval {

public:

  JetRecoEval(PHCompositeNode *topNode,
	      std::string recojetname,
	      std::string truthjetname);
  virtual ~JetRecoEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}

  JetTruthEval*     get_truth_eval()         {return &_jettrutheval;}
  SvtxEvalStack*    get_svtx_eval_stack()    {return _jetrutheval.get_svtx_eval_stack();}
  CaloEvalStack*    get_cemc_eval_stack()    {return _jetrutheval.get_cemc_eval_stack();}
  CaloEvalStack*    get_hcalin_eval_stack()  {return _jetrutheval.get_hcalin_eval_stack();}
  CaloEvalStack*    get_hcalout_eval_stack() {return _jetrutheval.get_hcalout_eval_stack();}

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

  JetTruthEval _jettrutheval;
  std::string _recojetname;
  std::string _truthjetname;
    
  JetMap* _recojets;
  JetMap* _truthjets;
  
  bool                                    _do_cache;
  std::map<Jet*,std::set<PHG4Hit*> >      _cache_all_truth_hits;
  std::map<Jet*,std::set<Jet*> >          _cache_all_truth_jets;
  std::map<Jet*,Jet* >                    _cache_max_truth_jet_by_energy;
  std::map<Jet*,std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<Jet*,std::set<Jet*> >          _cache_all_jets_from;
  std::map<Jet*,Jet* >                    _cache_best_jet_from;
  std::map<std::pair<Jet*,Jet*>,float>    _cache_get_energy_contribution;
};

#endif // __SVTXHITEVAL_H__
