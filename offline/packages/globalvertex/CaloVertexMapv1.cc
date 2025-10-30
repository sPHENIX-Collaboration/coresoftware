#include "CaloVertexMapv1.h"

#include "CaloVertex.h"

#include <iterator>  // for reverse_iterator
#include <utility>   // for pair, make_pair

CaloVertexMapv1::~CaloVertexMapv1()
{
  CaloVertexMapv1::clear();
}

void CaloVertexMapv1::identify(std::ostream& os) const
{
  os << "CaloVertexMapv1: size = " << _map.size() << std::endl;
  return;
}

void CaloVertexMapv1::clear()
{
  for (auto& iter : _map)
  {
    delete iter.second;
  }
  _map.clear();
  return;
}

const CaloVertex* CaloVertexMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

CaloVertex* CaloVertexMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

CaloVertex* CaloVertexMapv1::insert(CaloVertex* vertex)
{
  unsigned int index = 0;
  if (!_map.empty())
  {
    index = _map.rbegin()->first + 1;
  }
  _map.insert(std::make_pair(index, vertex));
  _map[index]->set_id(index);
  return _map[index];
}
