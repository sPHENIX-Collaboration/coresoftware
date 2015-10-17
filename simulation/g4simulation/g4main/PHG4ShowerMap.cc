#include "PHG4ShowerMap.h"

#include "PHG4Shower.h"

using namespace std;

ClassImp(PHG4ShowerMap)

PHG4ShowerMap::PHG4ShowerMap()
: _map() {
}

PHG4ShowerMap::~PHG4ShowerMap() {
  clear();
}

void PHG4ShowerMap::identify(ostream& os) const {
  os << "PHG4ShowerMap: size = " << _map.size() << endl;
  return;  
}

void PHG4ShowerMap::clear() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    delete iter->second;
  }
  _map.clear();
  return;
}

const PHG4Shower* PHG4ShowerMap::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

PHG4Shower* PHG4ShowerMap::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

PHG4Shower* PHG4ShowerMap::insert(PHG4Shower* shower) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , shower ));
  _map[index]->set_id(index);
  return _map[index];
}
