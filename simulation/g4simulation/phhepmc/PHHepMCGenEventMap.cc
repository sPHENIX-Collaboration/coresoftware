#include "PHHepMCGenEventMap.h"

#include "PHHepMCGenEvent.h"

using namespace std;

ClassImp(PHHepMCGenEventMap)

PHHepMCGenEventMap::PHHepMCGenEventMap()
: _map() {
}

PHHepMCGenEventMap::PHHepMCGenEventMap(const PHHepMCGenEventMap& eventmap)
  : _map() {  
  for (ConstIter iter = eventmap.begin();
       iter != eventmap.end();
       ++iter) {
    const PHHepMCGenEvent *event = iter->second;
    _map.insert(make_pair(event->get_embedding_id(),event->Clone()));
  }  
}

PHHepMCGenEventMap& PHHepMCGenEventMap::operator=(const PHHepMCGenEventMap& eventmap) {
  Reset();
  for (ConstIter iter = eventmap.begin();
       iter != eventmap.end();
       ++iter) {
    const PHHepMCGenEvent *event = iter->second;
    _map.insert(make_pair(event->get_embedding_id(),event->Clone()));
  }  
  return *this;
}

PHHepMCGenEventMap::~PHHepMCGenEventMap() {
  Reset();
}

void PHHepMCGenEventMap::Reset() {
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter) {
    PHHepMCGenEvent *event = iter->second;
    delete event;
  }
  _map.clear();
}

void PHHepMCGenEventMap::identify(ostream& os) const {
  os << "PHHepMCGenEventMap: size = " << _map.size() << endl;
  return;  
}

const PHHepMCGenEvent* PHHepMCGenEventMap::get(int id) const {
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return NULL;  
  return iter->second;
}

PHHepMCGenEvent* PHHepMCGenEventMap::get(int id) {
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

PHHepMCGenEvent* PHHepMCGenEventMap::insert_active_event(const PHHepMCGenEvent* event) {
  unsigned int index = 1;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair( index , event->Clone() ));
  _map[index]->set_embedding_id(index);
  return _map[index];
}

PHHepMCGenEvent* PHHepMCGenEventMap::insert_background_event(const PHHepMCGenEvent* event) {
  unsigned int index = 0;
  if (!_map.empty()) index = _map.begin()->first - 1;
  _map.insert(make_pair( index , event->Clone() ));
  _map[index]->set_embedding_id(index);
  return _map[index];
}
