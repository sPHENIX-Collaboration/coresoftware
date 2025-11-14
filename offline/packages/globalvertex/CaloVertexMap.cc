#include "CaloVertexMap.h"

#include <map>

class CaloVertex;

namespace
{
  std::map<unsigned int, CaloVertex *> DummyCaloVertexMap;
}

CaloVertexMap::ConstIter CaloVertexMap::begin() const
{
  return DummyCaloVertexMap.end();
}

CaloVertexMap::ConstIter CaloVertexMap::find(unsigned int /*idkey*/) const
{
  return DummyCaloVertexMap.end();
}

CaloVertexMap::ConstIter CaloVertexMap::end() const
{
  return DummyCaloVertexMap.end();
}

CaloVertexMap::Iter CaloVertexMap::begin()
{
  return DummyCaloVertexMap.end();
}

CaloVertexMap::Iter CaloVertexMap::find(unsigned int /*idkey*/)
{
  return DummyCaloVertexMap.end();
}

CaloVertexMap::Iter CaloVertexMap::end()
{
  return DummyCaloVertexMap.end();
}
