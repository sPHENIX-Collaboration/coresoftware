#include "SvtxClusterMap.h"

#include "SvtxCluster.h"

using namespace std;

ClassImp(SvtxClusterMap)

SvtxClusterMap::SvtxClusterMap()
: _map() {
}

SvtxClusterMap::SvtxClusterMap(const SvtxClusterMap& clustermap)
  : _map() {  
  for (ConstIter iter = clustermap.begin();
       iter != clustermap.end();
       ++iter) {
    const SvtxCluster *cluster = iter->second;
    _map.insert(make_pair(cluster->get_id(),cluster->Clone()));
  }  
}

SvtxClusterMap& SvtxClusterMap::operator=(const SvtxClusterMap& clustermap) {
  Reset();
  for (ConstIter iter = clustermap.begin();
       iter != clustermap.end();
       ++iter) {
    const SvtxCluster *cluster = iter->second;
    _map.insert(make_pair(cluster->get_id(),cluster->Clone()));
  }  
  return *this;
}

SvtxClusterMap::~SvtxClusterMap() {
  Reset();
}

void SvtxClusterMap::Reset() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    SvtxCluster *cluster = iter->second;
    delete cluster;
  }
  _map.clear();
}

void SvtxClusterMap::identify(ostream& os) const {
  os << "SvtxClusterMap: size = " << _map.size() << endl;
  return;  
}

const SvtxCluster* SvtxClusterMap::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

SvtxCluster* SvtxClusterMap::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxCluster* SvtxClusterMap::insert(const SvtxCluster* clus) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , clus->Clone() ));
  _map[index]->set_id(index);
  return _map[index];
}
