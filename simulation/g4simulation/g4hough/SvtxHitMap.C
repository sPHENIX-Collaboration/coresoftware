#include "SvtxHitMap.h"

#include "SvtxHit.h"

#include <map>

using namespace std;

ClassImp(SvtxHitMap)

SvtxHitMap::SvtxHitMap()
: _map() {
}

SvtxHitMap::~SvtxHitMap() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    SvtxHit *hit = iter->second;
    delete hit;
  }  
  _map.clear();
}

void SvtxHitMap::identify(ostream& os) const {
  os << "SvtxHitMap: size = " << _map.size() << endl;
  return;  
}

const SvtxHit* SvtxHitMap::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

SvtxHit* SvtxHitMap::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxHit* SvtxHitMap::insert(const SvtxHit *hit) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , hit->Clone() ));
  _map[index]->set_id(index);
  return _map[index];
}
