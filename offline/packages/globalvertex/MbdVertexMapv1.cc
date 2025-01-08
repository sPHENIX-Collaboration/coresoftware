#include "MbdVertexMapv1.h"

#include "MbdVertex.h"
#include "MbdVertexMap.h"

#include <iterator>  // for reverse_iterator
#include <utility>   // for pair, make_pair

using namespace std;

MbdVertexMapv1::MbdVertexMapv1()
  : _map()
{
}

MbdVertexMapv1::~MbdVertexMapv1()
{
  MbdVertexMapv1::clear();
}

void MbdVertexMapv1::identify(ostream& os) const
{
  os << "MbdVertexMapv1: size = " << _map.size() << endl;
  return;
}

void MbdVertexMapv1::clear()
{
  for (auto& iter : _map)
  {
    delete iter.second;
  }
  _map.clear();
  return;
}

const MbdVertex* MbdVertexMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

MbdVertex* MbdVertexMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

MbdVertex* MbdVertexMapv1::insert(MbdVertex* clus)
{
  unsigned int index = 0;
  if (!_map.empty())
  {
    index = _map.rbegin()->first + 1;
  }
  _map.insert(make_pair(index, clus));
  _map[index]->set_id(index);
  return _map[index];
}
