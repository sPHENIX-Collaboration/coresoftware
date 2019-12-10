#ifndef G4EVAL_CALOEVALSTACK_H
#define G4EVAL_CALOEVALSTACK_H

#include "CaloRawClusterEval.h"
#include "CaloTruthEval.h"

#include <string>

class CaloRawTowerEval;
class PHCompositeNode;

// This user class provides pointers to the
// full set of calorimeter evaluators and
// protects the user from future introduction
// of new eval heirachies (new eval objects can
// be introduced without rewrites)

class CaloEvalStack
{
 public:
  CaloEvalStack(PHCompositeNode* topNode, const std::string& caloname);
  virtual ~CaloEvalStack() {}

  int get_caloid() { return get_truth_eval()->get_caloid(); }
  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache) { _clustereval.do_caching(do_cache); }
  void set_strict(bool strict) { _clustereval.set_strict(strict); }
  void set_verbosity(int verbosity) { _clustereval.set_verbosity(verbosity); }

  CaloRawClusterEval* get_rawcluster_eval() { return &_clustereval; }
  CaloRawTowerEval* get_rawtower_eval() { return _clustereval.get_rawtower_eval(); }
  CaloTruthEval* get_truth_eval() { return _clustereval.get_truth_eval(); }

  unsigned int get_errors() { return _clustereval.get_errors(); }

 private:
  CaloRawClusterEval _clustereval;  // right now this is the top-level eval, other evals nest underneath
};

#endif  // G4EVAL_CALOEVALSTACK_H
