#include "GlobalVertexMapv1.h"

#include "GlobalVertex.h"
#include "GlobalVertexMap.h"

#include <iterator>           // for reverse_iterator
#include <utility>            // for pair, make_pair

using namespace std;

GlobalVertexMapv1::GlobalVertexMapv1()
  : _map()
{
}

GlobalVertexMapv1::~GlobalVertexMapv1()
{
  clear();
}

void GlobalVertexMapv1::identify(ostream& os) const
{
  os << "GlobalVertexMapv1: size = " << _map.size() << endl;
  return;
}

void GlobalVertexMapv1::clear()
{
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter)
  {
    delete iter->second;
  }
  _map.clear();
  return;
}

const GlobalVertex* GlobalVertexMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

GlobalVertex* GlobalVertexMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

GlobalVertex* GlobalVertexMapv1::insert(GlobalVertex* clus)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair(index, clus));
  _map[index]->set_id(index);
  return _map[index];
}
