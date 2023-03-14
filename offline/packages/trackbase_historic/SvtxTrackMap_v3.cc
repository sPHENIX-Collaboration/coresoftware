#include "SvtxTrackMap_v3.h"

#include "SvtxTrack.h"

#include <phool/PHObject.h>  // for PHObject

#include <iterator>     // for reverse_iterator
#include <map>          // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>      // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>      // for pair, make_pair

SvtxTrackMap_v3::SvtxTrackMap_v3()
  : _map()
{
}

SvtxTrackMap_v3::SvtxTrackMap_v3(const SvtxTrackMap_v3& trackmap)
  : _map()
{
  for (ConstIter iter = trackmap.begin();
       iter != trackmap.end();
       ++iter)
  {
    auto track = static_cast<SvtxTrack*> (iter->second->CloneMe());
    _map.insert(std::make_pair(track->get_id(), track));
  }
}

SvtxTrackMap_v3& SvtxTrackMap_v3::operator=(const SvtxTrackMap_v3& trackmap)
{
  
  // do nothing if same  copying map onto itself
  if( &trackmap == this ) return *this;
  
  Reset();
  for (ConstIter iter = trackmap.begin();
       iter != trackmap.end();
       ++iter)
  {
    auto track = static_cast<SvtxTrack*> (iter->second->CloneMe());
    _map.insert(std::make_pair(track->get_id(), track));
  }
  return *this;
}

SvtxTrackMap_v3::~SvtxTrackMap_v3()
{
  Reset();
}

void SvtxTrackMap_v3::Reset()
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

void SvtxTrackMap_v3::identify(std::ostream& os) const
{
  os << "SvtxTrackMap_v3: size = " << _map.size() << std::endl;
  return;
}

const SvtxTrack* SvtxTrackMap_v3::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

SvtxTrack* SvtxTrackMap_v3::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

SvtxTrack* SvtxTrackMap_v3::insert(const SvtxTrack* track)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  auto copy = static_cast<SvtxTrack*>( track->CloneMe() );
  copy->set_id(index);

  const auto result = _map.insert(std::make_pair(index, copy ));
  if( !result.second ) 
  {
    std::cout << "SvtxTrackMap_v3::insert - duplicated key. track not inserted" << std::endl; 
    delete copy;
    return nullptr;
  } else {
    return copy;
  }
}

SvtxTrack* SvtxTrackMap_v3::insertWithKey(const SvtxTrack* track, unsigned int index)
{
  auto copy = static_cast<SvtxTrack*>( track->CloneMe() );
  copy->set_id(index);
  const auto result = _map.insert(std::make_pair(index, copy ));
  if( !result.second ) 
  {
    std::cout << "SvtxTrackMap_v3::insertWithKey - duplicated key. track not inserted" << std::endl; 
    delete copy;
    return nullptr;
  } else {
    return copy;
  }
}
