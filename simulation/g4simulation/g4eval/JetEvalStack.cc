#include "JetEvalStack.h"

JetEvalStack::JetEvalStack(PHCompositeNode* topNode,
                           const std::string& recojetname,
                           const std::string& truthjetname)
  : _recoeval(topNode, recojetname, truthjetname)
{
}

void JetEvalStack::next_event(PHCompositeNode* topNode)
{
  _recoeval.next_event(topNode);
}
