#include "SvtxVertexMap.h"

SvtxVertexMap::VertexMap DummyVertexMap;

SvtxVertexMap::ConstIter SvtxVertexMap::begin() const
{
  return DummyVertexMap.end();
}

SvtxVertexMap::ConstIter SvtxVertexMap::find(unsigned int /*idkey*/) const
{
  return DummyVertexMap.end();
}

SvtxVertexMap::ConstIter SvtxVertexMap::end() const
{
  return DummyVertexMap.end();
}


SvtxVertexMap::Iter SvtxVertexMap::begin()
{
  return DummyVertexMap.end();
}

SvtxVertexMap::Iter SvtxVertexMap::find(unsigned int /*idkey*/)
{
  return DummyVertexMap.end();
}

SvtxVertexMap::Iter SvtxVertexMap::end()
{
  return DummyVertexMap.end();
}

