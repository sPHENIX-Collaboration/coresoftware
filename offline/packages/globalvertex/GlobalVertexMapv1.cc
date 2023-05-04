#include "GlobalVertexMapv1.h"

#include "GlobalVertex.h"
#include "GlobalVertexMap.h"

#include <iterator>  // for reverse_iterator
#include <utility>   // for pair, make_pair

GlobalVertexMapv1::~GlobalVertexMapv1()
{
  clear();
}

void GlobalVertexMapv1::identify(std::ostream& os) const
{
  os << "GlobalVertexMapv1: size = " << _map.size() << std::endl;
  for (auto& m : _map)
  {
    m.second->identify(os);
  }
  return;
}

void GlobalVertexMapv1::clear()
{
  for (const auto& iter : _map)
  {
    delete iter.second;
  }
  _map.clear();
  return;
}

const GlobalVertex* GlobalVertexMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

GlobalVertex* GlobalVertexMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

GlobalVertex* GlobalVertexMapv1::insert(GlobalVertex* clus)
{
  unsigned int index = clus->get_id();
  auto ret = _map.insert(std::make_pair(index, clus));
  return ret.first->second;
}
