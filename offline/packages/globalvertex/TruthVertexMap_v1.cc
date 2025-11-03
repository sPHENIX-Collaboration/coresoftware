#include "TruthVertexMap_v1.h"

TruthVertexMap_v1::~TruthVertexMap_v1()
{
  clear();
}

void TruthVertexMap_v1::identify(std::ostream& os) const
{
  os << "TruthVertexMap_v1 with " << _map.size() << " entries" << std::endl;
}

void TruthVertexMap_v1::clear()
{
  for (auto& [id, vtx] : _map)
  {
    delete vtx;
  }
  _map.clear();
}

const TruthVertex* TruthVertexMap_v1::get(unsigned int idkey) const
{
  auto it = _map.find(idkey);
  return (it != _map.end()) ? it->second : nullptr;
}

TruthVertex* TruthVertexMap_v1::get(unsigned int idkey)
{
  auto it = _map.find(idkey);
  return (it != _map.end()) ? it->second : nullptr;
}

TruthVertex* TruthVertexMap_v1::insert(TruthVertex* vertex)
{
  if (!vertex) return nullptr;

  unsigned int id = vertex->get_id();
  auto [it, inserted] = _map.insert({id, vertex});

  if (!inserted)
  {
    delete it->second;
    it->second = vertex;
  }

  return vertex;
}

size_t TruthVertexMap_v1::erase(unsigned int idkey)
{
  auto it = _map.find(idkey);
  if (it != _map.end())
  {
    delete it->second;
    _map.erase(it);
    return 1;
  }
  return 0;
}
