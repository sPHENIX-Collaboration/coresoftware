#include "InttVertexMap.h"

// class BbcVertex;

std::map<unsigned int, InttVertex*> DummyInttVertexMap;

InttVertexMap::ConstIter InttVertexMap::begin() const
{
  return DummyInttVertexMap.end();
}

InttVertexMap::ConstIter InttVertexMap::find(unsigned int /*idkey*/) const
{
  return DummyInttVertexMap.end();
}

InttVertexMap::ConstIter InttVertexMap::end() const
{
  return DummyInttVertexMap.end();
}

InttVertexMap::Iter InttVertexMap::begin()
{
  return DummyInttVertexMap.end();
}

InttVertexMap::Iter InttVertexMap::find(unsigned int /*idkey*/)
{
  return DummyInttVertexMap.end();
}

InttVertexMap::Iter InttVertexMap::end()
{
  return DummyInttVertexMap.end();
}
