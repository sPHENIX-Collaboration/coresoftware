#include "BbcVertexMap.h"

class BbcVertex;

std::map<unsigned int, BbcVertex*> DummyBbcVertexMap;

BbcVertexMap::ConstIter BbcVertexMap::begin() const
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::ConstIter BbcVertexMap::find(unsigned int /*idkey*/) const
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::ConstIter BbcVertexMap::end() const
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::Iter BbcVertexMap::begin()
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::Iter BbcVertexMap::find(unsigned int /*idkey*/)
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::Iter BbcVertexMap::end()
{
  return DummyBbcVertexMap.end();
}
