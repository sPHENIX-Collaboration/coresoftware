#include "GlobalVertexMap_v1.h"

#include "GlobalVertexMap.h"
#include "GlobalVertex.h"

using namespace std;

ClassImp(GlobalVertexMap_v1)

GlobalVertexMap_v1::GlobalVertexMap_v1()
: _map() {
}

GlobalVertexMap_v1::~GlobalVertexMap_v1() {
  clear();
}

void GlobalVertexMap_v1::identify(ostream& os) const {
  os << "GlobalVertexMap_v1: size = " << _map.size() << endl;
  return;  
}

void GlobalVertexMap_v1::clear() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    delete iter->second;
  }
  _map.clear();
  return;
}

const GlobalVertex* GlobalVertexMap_v1::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

GlobalVertex* GlobalVertexMap_v1::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

GlobalVertex* GlobalVertexMap_v1::insert(GlobalVertex* clus) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , clus ));
  _map[index]->set_id(index);
  return _map[index];
}
