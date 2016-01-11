
#ifndef __JETEVALSTACK_H__
#define __JETEVALSTACK_H__

#include "JetRecoEval.h"
#include "JetTruthEval.h"

#include "SvtxEvalStack.h"
#include "CaloEvalStack.h"

#include <phool/PHCompositeNode.h>

// This user class provides pointers to the
// full set of jet evaluators and
// protects the user from future introduction
// of new eval heirachies (new eval objects can
// be introduced without rewrites)

class JetEvalStack {

public:

  JetEvalStack(PHCompositeNode *topNode,
	       std::string recojetname,
	       std::string truthjetname);
  virtual ~JetEvalStack() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_recoeval.do_caching(do_cache);}
  void set_strict(bool strict) {_recoeval.set_strict(strict);}
  void set_verbosity(int verbosity) {_recoeval.set_verbosity(verbosity);}
  
  JetRecoEval*   get_reco_eval() {return &_recoeval;}
  JetTruthEval*  get_truth_eval() {return _recoeval.get_truth_eval();}

  SvtxEvalStack* get_stvx_eval_stack()    {return _recoeval.get_svtx_eval_stack();}
  CaloEvalStack* get_cemc_eval_stack()    {return _recoeval.get_cemc_eval_stack();}
  CaloEvalStack* get_hcalin_eval_stack()  {return _recoeval.get_hcalin_eval_stack();}
  CaloEvalStack* get_hcalout_eval_stack() {return _recoeval.get_hcalout_eval_stack();}
  CaloEvalStack* get_femc_eval_stack()    {return _recoeval.get_femc_eval_stack();}
  CaloEvalStack* get_fhcal_eval_stack()   {return _recoeval.get_fhcal_eval_stack();}
  
private:
  JetRecoEval _recoeval; 
};

#endif // __JETEVALSTACK_H__
