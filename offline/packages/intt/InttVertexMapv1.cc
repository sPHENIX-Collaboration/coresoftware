#include "InttVertexMapv1.h"

#include "InttVertex.h"

#include <iterator>  // for reverse_iterator
#include <utility>   // for pair, make_pair

InttVertexMapv1::~InttVertexMapv1()
{
  InttVertexMapv1::clear();
}

void InttVertexMapv1::identify(std::ostream& os) const
{
  os << "InttVertexMapv1: size = " << _map.size() << std::endl;
  return;
}

void InttVertexMapv1::clear()
{
  for (auto& iter : _map)
  {
    delete iter.second;
  }
  _map.clear();
  return;
}

const InttVertex* InttVertexMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

InttVertex* InttVertexMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

InttVertex* InttVertexMapv1::insert(InttVertex* clus)
{
  unsigned int index = 0;
  if (!_map.empty())
  {
    index = _map.rbegin()->first + 1;
  }
  _map.insert(std::make_pair(index, clus));
  _map[index]->set_id(index);
  return _map[index];
}
