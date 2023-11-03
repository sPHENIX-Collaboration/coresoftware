#include "MbdVertexMap.h"

class MbdVertex;

std::map<unsigned int, MbdVertex*> DummyMbdVertexMap;

MbdVertexMap::ConstIter MbdVertexMap::begin() const
{
  return DummyMbdVertexMap.end();
}

MbdVertexMap::ConstIter MbdVertexMap::find(unsigned int /*idkey*/) const
{
  return DummyMbdVertexMap.end();
}

MbdVertexMap::ConstIter MbdVertexMap::end() const
{
  return DummyMbdVertexMap.end();
}

MbdVertexMap::Iter MbdVertexMap::begin()
{
  return DummyMbdVertexMap.end();
}

MbdVertexMap::Iter MbdVertexMap::find(unsigned int /*idkey*/)
{
  return DummyMbdVertexMap.end();
}

MbdVertexMap::Iter MbdVertexMap::end()
{
  return DummyMbdVertexMap.end();
}
