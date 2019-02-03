#include "CaloEvalStack.h"

CaloEvalStack::CaloEvalStack(PHCompositeNode* topNode, const std::string& caloname)
  : _clustereval(topNode, caloname)
{
}

void CaloEvalStack::next_event(PHCompositeNode* topNode)
{
  _clustereval.next_event(topNode);
}
