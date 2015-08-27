
#include "SvtxEvalStack.h"

#include "SvtxVertexEval.h"
#include "SvtxTruthEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>

#include <string>

using namespace std;

SvtxEvalStack::SvtxEvalStack(PHCompositeNode* topNode)
  : _vertexeval(topNode),
    _trutheval(topNode) {
}

void SvtxEvalStack::next_event(PHCompositeNode* topNode) {
  _vertexeval.next_event(topNode);
  _trutheval.next_event(topNode);
}

