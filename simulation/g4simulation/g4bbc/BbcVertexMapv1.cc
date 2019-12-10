#include "BbcVertexMapv1.h"

#include "BbcVertex.h"
#include "BbcVertexMap.h"

#include <iterator>        // for reverse_iterator
#include <utility>         // for pair, make_pair

using namespace std;

BbcVertexMapv1::BbcVertexMapv1()
  : _map()
{
}

BbcVertexMapv1::~BbcVertexMapv1()
{
  clear();
}

void BbcVertexMapv1::identify(ostream& os) const
{
  os << "BbcVertexMapv1: size = " << _map.size() << endl;
  return;
}

void BbcVertexMapv1::clear()
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

const BbcVertex* BbcVertexMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

BbcVertex* BbcVertexMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

BbcVertex* BbcVertexMapv1::insert(BbcVertex* clus)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair(index, clus));
  _map[index]->set_id(index);
  return _map[index];
}
