
#ifndef __CALOEVALSTACK_H__
#define __CALOEVALSTACK_H__

#include "CaloRawClusterEval.h"
#include "CaloRawTowerEval.h"
#include "CaloTruthEval.h"

#include <phool/PHCompositeNode.h>

#include <string>

// This user class provides pointers to the
// full set of calorimeter evaluators and
// protects the user from future introduction
// of new eval heirachies (new eval objects can
// be introduced without rewrites)

class CaloEvalStack {

public:

  CaloEvalStack(PHCompositeNode *topNode, std::string caloname);
  virtual ~CaloEvalStack() {}

  void next_event(PHCompositeNode *topNode);

  CaloRawClusterEval* get_rawcluster_eval() {return &_clustereval;}
  CaloRawTowerEval*   get_rawtower_eval() {return _clustereval->get_rawtower_eval();}
  CaloTruthEval*      get_truth_eval() {return _clustereval->get_truth_eval();}
  
private:
  CaloRawClusterEval _clustereval; // right now this is the top-level eval, other evals nest underneath
};

#endif // __CALOEVALSTACK_H__
