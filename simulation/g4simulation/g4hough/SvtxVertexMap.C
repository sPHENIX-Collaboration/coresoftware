#include "SvtxVertexMap.h"

#include "SvtxVertex.h"

using namespace std;

ClassImp(SvtxVertexMap)

SvtxVertexMap::SvtxVertexMap()
: _map() {
}

SvtxVertexMap::~SvtxVertexMap() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    SvtxVertex *vertex = iter->second;
    delete vertex;
  }  
  _map.clear();
}

void SvtxVertexMap::identify(ostream& os) const {
  os << "SvtxVertexMap: size = " << _map.size() << endl;
  return;  
}

const SvtxVertex* SvtxVertexMap::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

SvtxVertex* SvtxVertexMap::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxVertex* SvtxVertexMap::insert(const SvtxVertex &clus) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , new SvtxVertex(clus) ));
  _map[index]->set_id(index);
  return _map[index];
}
