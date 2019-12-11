#ifndef G4EVAL_JETEVALSTACK_H
#define G4EVAL_JETEVALSTACK_H

#include "JetRecoEval.h"

#include <string>

class CaloEvalStack;
class JetTruthEval;
class SvtxEvalStack;
class PHCompositeNode;

// This user class provides pointers to the
// full set of jet evaluators and
// protects the user from future introduction
// of new eval heirachies (new eval objects can
// be introduced without rewrites)

class JetEvalStack
{
 public:
  JetEvalStack(PHCompositeNode* topNode,
               const std::string& recojetname,
               const std::string& truthjetname);
  virtual ~JetEvalStack() {}

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache) { _recoeval.do_caching(do_cache); }
  void set_strict(bool strict) { _recoeval.set_strict(strict); }
  void set_verbosity(int verbosity) { _recoeval.set_verbosity(verbosity); }

  JetRecoEval* get_reco_eval() { return &_recoeval; }
  JetTruthEval* get_truth_eval() { return _recoeval.get_truth_eval(); }

  SvtxEvalStack* get_stvx_eval_stack() { return _recoeval.get_svtx_eval_stack(); }
  CaloEvalStack* get_cemc_eval_stack() { return _recoeval.get_cemc_eval_stack(); }
  CaloEvalStack* get_hcalin_eval_stack() { return _recoeval.get_hcalin_eval_stack(); }
  CaloEvalStack* get_hcalout_eval_stack() { return _recoeval.get_hcalout_eval_stack(); }
  CaloEvalStack* get_femc_eval_stack() { return _recoeval.get_femc_eval_stack(); }
  CaloEvalStack* get_fhcal_eval_stack() { return _recoeval.get_fhcal_eval_stack(); }

 private:
  JetRecoEval _recoeval;
};

#endif  // G4EVAL_JETEVALSTACK_H
