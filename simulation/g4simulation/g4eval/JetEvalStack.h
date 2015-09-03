
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

  JetEvalStack(PHCompositeNode *topNode, std::string jetnode);
  virtual ~JetEvalStack() {}

  void next_event(PHCompositeNode *topNode);

  JetRecoEval*   get_reco_eval() {return &_recoeval;}
  JetTruthEval*  get_truth_eval() {return _recoeval->get_truth_eval;}

  SvtxEvalStack* get_stvx_eval_stack()                     {return _recoeval->get_svtx_eval_stack();}
  CaloEvalStack* get_calo_eval_stack(std::string caloname) {return _recoeval->get_calo_eval_stack(caloname);}
  
private:
  JetRecoEval _recoeval; 
};

#endif // __JETEVALSTACK_H__
