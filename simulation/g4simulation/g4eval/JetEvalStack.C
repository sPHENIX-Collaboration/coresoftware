
#include "JetEvalStack.h"

#include "JetRecoEval.h"

#include <phool/PHCompositeNode.h>

#include <string>

using namespace std;

JetEvalStack::JetEvalStack(PHCompositeNode* topNode,
			   std::string recojetname,
			   std::string truthjetname)
  : _recoeval(topNode,recojetname,truthjetname) {}

void JetEvalStack::next_event(PHCompositeNode* topNode) {
  _recoeval.next_event(topNode);
}

