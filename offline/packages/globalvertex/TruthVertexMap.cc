#include "TruthVertexMap.h"

class TruthVertex;

TruthVertexMap::ConstIter TruthVertexMap::begin() const { return ConstIter(); }
TruthVertexMap::ConstIter TruthVertexMap::find(unsigned int) const { return ConstIter(); }
TruthVertexMap::ConstIter TruthVertexMap::end() const { return ConstIter(); }

TruthVertexMap::Iter TruthVertexMap::begin() { return Iter(); }
TruthVertexMap::Iter TruthVertexMap::find(unsigned int) { return Iter(); }
TruthVertexMap::Iter TruthVertexMap::end() { return Iter(); }