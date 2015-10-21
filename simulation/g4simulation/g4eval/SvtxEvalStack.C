
#include "SvtxEvalStack.h"

#include "SvtxVertexEval.h"

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <string>

using namespace std;

SvtxEvalStack::SvtxEvalStack(PHCompositeNode* topNode)
  : _vertexeval(topNode) {
}

void SvtxEvalStack::next_event(PHCompositeNode* topNode) {
  _vertexeval.next_event(topNode);
}

