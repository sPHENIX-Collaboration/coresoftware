#include "BbcVertexMap.h"

#include "BbcVertex.h"

using namespace std;

ClassImp(BbcVertexMap)

BbcVertexMap::BbcVertexMap()
: _map() {
}

BbcVertexMap::~BbcVertexMap() {
  clear();
}

void BbcVertexMap::identify(ostream& os) const {
  os << "BbcVertexMap: size = " << _map.size() << endl;
  return;  
}

void BbcVertexMap::clear() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    delete iter->second;
  }
  _map.clear();
  return;
}

const BbcVertex* BbcVertexMap::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

BbcVertex* BbcVertexMap::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

BbcVertex* BbcVertexMap::insert(BbcVertex* clus) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , clus ));
  _map[index]->set_id(index);
  return _map[index];
}
