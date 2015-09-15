#include "SvtxTrackMap.h"

#include "SvtxTrack.h"

using namespace std;

ClassImp(SvtxTrackMap)

SvtxTrackMap::SvtxTrackMap()
: _map() {
}

SvtxTrackMap::~SvtxTrackMap() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    SvtxTrack *track = iter->second;
    delete track;
  }  
  _map.clear();
}

void SvtxTrackMap::identify(ostream& os) const {
  os << "SvtxTrackMap: size = " << _map.size() << endl;
  return;  
}

const SvtxTrack* SvtxTrackMap::get(unsigned int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

SvtxTrack* SvtxTrackMap::get(unsigned int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxTrack* SvtxTrackMap::insert(const SvtxTrack &track) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , new SvtxTrack(track) ));
  _map[index]->set_id(index);
  return _map[index];
}
