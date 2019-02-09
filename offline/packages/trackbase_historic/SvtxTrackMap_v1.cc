#include "SvtxTrackMap_v1.h"

#include "SvtxTrack.h"

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
    const SvtxTrack* track = iter->second;
    _map.insert(make_pair(track->get_id(), track->Clone()));
  }
}

SvtxTrackMap_v1& SvtxTrackMap_v1::operator=(const SvtxTrackMap_v1& trackmap)
{
  Reset();
  for (ConstIter iter = trackmap.begin();
       iter != trackmap.end();
       ++iter)
  {
    const SvtxTrack* track = iter->second;
    _map.insert(make_pair(track->get_id(), track->Clone()));
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
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxTrack* SvtxTrackMap_v1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return NULL;
  return iter->second;
}

SvtxTrack* SvtxTrackMap_v1::insert(const SvtxTrack* track)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(make_pair(index, track->Clone()));
  _map[index]->set_id(index);
  return _map[index];
}
