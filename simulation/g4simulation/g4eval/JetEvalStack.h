
#ifndef __JETEVALSTACK_H__
#define __JETEVALSTACK_H__

#include "JetRecoEval.h"
#include "JetTruthEval.h"

#include <phool/PHCompositeNode.h>

// This user class provides pointers to the
// full set of jet evaluators and
// protects the user from future introduction
// of new eval heirachies (new eval objects can
// be introduced without rewrites)

class JetEvalStack {

public:

  JetEvalStack(PHCompositeNode *topNode);
  virtual ~JetEvalStack() {}

  void next_event(PHCompositeNode *topNode);

  JetRecoEval*  get_reco_eval() {return &_recoeval;}
  JetTruthEval* get_truth_eval() {return &_trutheval;}
  
private:
  JetRecoEval _recoeval; // right now this is the top-level eval, other evals nest underneath
  JetTruthEval  _trutheval;  // except this one
};

#endif // __JETEVALSTACK_H__
