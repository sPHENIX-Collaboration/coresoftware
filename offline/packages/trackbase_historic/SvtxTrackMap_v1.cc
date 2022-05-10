#include "SvtxTrackMap_v1.h"

#include "SvtxTrack.h"

#include <phool/PHObject.h>  // for PHObject

#include <iterator>     // for reverse_iterator
#include <map>          // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>      // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>      // for pair, make_pair

using namespace std;


SvtxTrackMap_v1::SvtxTrackMap_v1()
  : _map()
{
}

SvtxTrackMap_v1::SvtxTrackMap_v1(const SvtxTrackMap_v1& trackmap)
  : _map()
{
  for (ConstIter iter = trackmap.begin();
       iter != trackmap.end();
       ++iter)
  {
    SvtxTrack* track = dynamic_cast<SvtxTrack*> (iter->second->CloneMe());
    _map.insert(make_pair(track->get_id(), track));
  }
}

SvtxTrackMap_v1& SvtxTrackMap_v1::operator=(const SvtxTrackMap_v1& trackmap)
{
  Reset();
  for (ConstIter iter = trackmap.begin();
       iter != trackmap.end();
       ++iter)
  {
    SvtxTrack* track = dynamic_cast<SvtxTrack*> (iter->second->CloneMe());
    _map.insert(make_pair(track->get_id(), track));
  }
  return *this;
}

SvtxTrackMap_v1::~SvtxTrackMap_v1()
{
  Reset();
}

void SvtxTrackMap_v1::Reset()
{
  for (Iter iter = _map.begin();
       iter != _map.end();
       ++iter)
  {
    SvtxTrack* track = iter->second;
    delete track;
  }
  _map.clear();
}

void SvtxTrackMap_v1::identify(ostream& os) const
{
  os << "SvtxTrackMap_v1: size = " << _map.size() << endl;
  return;
}

const SvtxTrack* SvtxTrackMap_v1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

SvtxTrack* SvtxTrackMap_v1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

SvtxTrack* SvtxTrackMap_v1::insert(const SvtxTrack* track)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair(index, dynamic_cast<SvtxTrack*> (track->CloneMe())));
  _map[index]->set_id(index);
  return _map[index];
}

size_t SvtxTrackMap_v1::erase(unsigned int idkey)
{
  const auto iter = _map.find(idkey);
  if( iter != _map.end() ) 
  { 
  
    delete iter->second;
    _map.erase( iter );
    return 1;
  
  } else return 0;
}

