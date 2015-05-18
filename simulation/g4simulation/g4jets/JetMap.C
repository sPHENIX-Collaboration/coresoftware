#include "JetMap.h"

#include "Jet.h"

#include <cmath>

using namespace std;

ClassImp(JetMap)

JetMap::JetMap()
: _algo(Jet::NONE),
  _par(NAN),
  _src(),
  _map() {
}

JetMap::~JetMap() {
  JetMap::Reset();
}

void JetMap::Reset() {
  _algo = Jet::NONE;  
  _par = NAN;
  _src.clear();
//  _map.clear();

  while(_map.begin() != _map.end())
    {
      delete _map.begin()->second;
      _map.erase(_map.begin());
    }
}

void JetMap::identify(ostream& os) const {
  os << "JetMap: size = " << _map.size() << endl;
  return;  
}

const Jet* JetMap::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

Jet* JetMap::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

Jet* JetMap::insert(Jet* jet) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , jet ));
  _map[index]->set_id(index);
  return (_map[index]);
}
