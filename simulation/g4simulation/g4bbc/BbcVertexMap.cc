#include "BbcVertexMap.h"

std::map<unsigned int, BbcVertex*> DummyBbcVertexMap;

BbcVertexMap::ConstIter BbcVertexMap::begin() const
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::ConstIter BbcVertexMap::find(unsigned int) const
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::ConstIter BbcVertexMap::end() const
{
  return DummyBbcVertexMap.end();
}


BbcVertexMap::Iter BbcVertexMap:: begin()
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::Iter BbcVertexMap::find(unsigned int)
{
  return DummyBbcVertexMap.end();
}

BbcVertexMap::Iter BbcVertexMap::end()
{
  return DummyBbcVertexMap.end();
}
