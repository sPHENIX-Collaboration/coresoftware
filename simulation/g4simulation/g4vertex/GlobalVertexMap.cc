#include "GlobalVertexMap.h"

std::map<unsigned int, GlobalVertex*> DummyGlobalVertexMap;

GlobalVertexMap::ConstIter GlobalVertexMap::begin() const
{
  return DummyGlobalVertexMap.end();
}

GlobalVertexMap::ConstIter GlobalVertexMap::find(unsigned int) const
{
  return DummyGlobalVertexMap.end();
}

GlobalVertexMap::ConstIter GlobalVertexMap::end() const
{
  return DummyGlobalVertexMap.end();
}

GlobalVertexMap::Iter GlobalVertexMap::begin()
{
  return DummyGlobalVertexMap.end();
}

GlobalVertexMap::Iter GlobalVertexMap::find(unsigned int)
{
  return DummyGlobalVertexMap.end();
}

GlobalVertexMap::Iter GlobalVertexMap::end()
{
  return DummyGlobalVertexMap.end();
}

