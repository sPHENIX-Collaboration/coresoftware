#include "BbcVertexMapv1.h"

#include "BbcVertex.h"
#include "BbcVertexMap.h"

#include <iterator>  // for reverse_iterator
#include <utility>   // for pair, make_pair

BbcVertexMapv1::~BbcVertexMapv1()
{
  BbcVertexMapv1::clear();
}

void BbcVertexMapv1::identify(std::ostream& os) const
{
  os << "BbcVertexMapv1: size = " << _map.size() << std::endl;
  return;
}

void BbcVertexMapv1::clear()
{
  for (auto& iter : _map)
  {
    delete iter.second;
  }
  _map.clear();
  return;
}

const BbcVertex* BbcVertexMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

BbcVertex* BbcVertexMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

BbcVertex* BbcVertexMapv1::insert(BbcVertex* vertex)
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
