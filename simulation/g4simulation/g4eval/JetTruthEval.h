
#ifndef __JETTRUTHEVAL_H__
#define __JETTRUTHEVAL_H__

#include "SvtxEvalStack.h"
#include "CaloEvalStack.h"

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4jets/Jet.h>

#include <map>
#include <string>

class JetTruthEval {

public:

  JetTruthEval(PHCompositeNode* topNode);
  virtual ~JetTruthEval() {}

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}

  SvtxEvalStack* get_svtx_eval_stack() {return &_svtxevalstack;}
  CaloEvalStack* get_cemc_eval_stack() {return &_cemcevalstack;}
  CaloEvalStack* get_hcalin_eval_stack() {return &_hcalinevalstack;}
  CaloEvalStack* get_hcalout_eval_stack() {return &_hcaloutevalstack;}
  
  std::set<PHG4Particle*> all_truth_particles(Jet* truthjet);
  std::set<PHG4Hit*>      all_truth_hits(Jet* truthjet);
  
private:

  void get_node_pointers(PHCompositeNode* topNode);

  SvtxEvalStack _svtxevalstack;
  CaloEvalStack _cemcevalstack;
  CaloEvalStack _hcalinevalstack;
  CaloEvalStack _hcaloutevalstack;

  PHG4TruthInfoContainer* _truthinfo;
  
  bool                                    _do_cache;
  std::map<Jet*,std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<Jet*,std::set<PHG4Hit*> >      _cache_all_truth_hits;
};

#endif // __JETTRUTHEVAL_H__
