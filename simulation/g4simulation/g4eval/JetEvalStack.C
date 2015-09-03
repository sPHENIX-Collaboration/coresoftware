
#include "JetEvalStack.h"

#include "JetRecoEval.h"
#include "JetTruthEval.h"

#include "SvtxEvalStack.h"
#include "CaloEvalStack.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>

#include <string>

using namespace std;

JetEvalStack::JetEvalStack(PHCompositeNode* topNode, std::string jetnode)
  : _recoeval(topNode,jetnode) {}

void JetEvalStack::next_event(PHCompositeNode* topNode) {
  _recoeval.next_event(topNode);
}

