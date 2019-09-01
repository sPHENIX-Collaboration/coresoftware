#include "SvtxHitMap_v1.h"

#include "SvtxHit.h"

#include <cassert>
#include <iterator>   // for reverse_iterator
#include <map>
#include <utility>    // for pair, make_pair

using namespace std;


    SvtxHitMap_v1::SvtxHitMap_v1()
  : _map()
{
}

SvtxHitMap_v1::SvtxHitMap_v1(const SvtxHitMap_v1& hitmap)
  : _map()
{
  for (ConstIter iter = hitmap.begin();
       iter != hitmap.end();
       ++iter)
  {
    SvtxHit* hit = dynamic_cast<SvtxHit *> (iter->second->CloneMe());
    _map.insert(make_pair(hit->get_id(), hit));
  }
}

SvtxHitMap_v1& SvtxHitMap_v1::operator=(const SvtxHitMap_v1& hitmap)
{
  Reset();
  for (ConstIter iter = hitmap.begin();
       iter != hitmap.end();
       ++iter)
  {
    SvtxHit* hit = dynamic_cast<SvtxHit *> (iter->second->CloneMe());
    _map.insert(make_pair(hit->get_id(), hit));
  }
  return *this;
}

SvtxHitMap_v1::~SvtxHitMap_v1()
{
  Reset();
}

void SvtxHitMap_v1::Reset()
{
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter)
  {
    SvtxHit* hit = iter->second;
    delete hit;
  }
  _map.clear();
}

void SvtxHitMap_v1::identify(ostream& os) const
{
  os << "SvtxHitMap_v1: size = " << _map.size() << endl;
  return;
}

const SvtxHit* SvtxHitMap_v1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxHit* SvtxHitMap_v1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxHit* SvtxHitMap_v1::insert(const SvtxHit* hit)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair(index, dynamic_cast<SvtxHit *> (hit->CloneMe())));
  _map[index]->set_id(index);
  return _map[index];
}
