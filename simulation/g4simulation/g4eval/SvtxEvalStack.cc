#include "SvtxEvalStack.h"

SvtxEvalStack::SvtxEvalStack(PHCompositeNode* topNode)
  : _vertexeval(topNode)
{
}

void SvtxEvalStack::next_event(PHCompositeNode* topNode)
{
  _vertexeval.next_event(topNode);
}
