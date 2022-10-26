#include "SvtxTrackCaloClusterMap_v1.h"

#include "SvtxTrack.h"

#include <phool/PHObject.h>  // for PHObject

#include <iterator>  // for reverse_iterator
#include <map>       // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>   // for pair, make_pair

using namespace std;

SvtxTrackCaloClusterMap_v1::SvtxTrackCaloClusterMap_v1()
  : _map()
{
}

SvtxTrackCaloClusterMap_v1::SvtxTrackCaloClusterMap_v1(const SvtxTrackCaloClusterMap_v1& map)
  : _map()
{
  for (ConstIter iter = map.begin();
       iter != map.end();
       ++iter)
  {
    std::vector<RawCluster*> clus = iter->second;
    _map.insert(make_pair(iter->first, clus));
  }
}

SvtxTrackCaloClusterMap_v1& SvtxTrackCaloClusterMap_v1::operator=(const SvtxTrackCaloClusterMap_v1& map)
{
  Reset();
  for (ConstIter iter = map.begin();
       iter != map.end();
       ++iter)
  {
    std::vector<RawCluster*> clus = iter->second;
    _map.insert(make_pair(iter->first, clus));
  }
  return *this;
}

SvtxTrackCaloClusterMap_v1::~SvtxTrackCaloClusterMap_v1()
{
  Reset();
}

void SvtxTrackCaloClusterMap_v1::Reset()
{
  /// Pointers to SvtxTrack and RawCluster will be
  /// deleted by their respective containers already
  _map.clear();
}

void SvtxTrackCaloClusterMap_v1::identify(ostream& os) const
{
  os << "SvtxTrackCaloClusterMap_v1: size = " << _map.size() << endl;
  return;
}

const std::vector<RawCluster*> SvtxTrackCaloClusterMap_v1::get(SvtxTrack* track) const
{
  std::vector<RawCluster*> dummy;
  ConstIter iter = _map.find(track);
  if (iter == _map.end())
  {
    return dummy;
  }
  return iter->second;
}

std::vector<RawCluster*> SvtxTrackCaloClusterMap_v1::get(SvtxTrack* track)
{
  std::vector<RawCluster*> dummy;
  Iter iter = _map.find(track);
  if (iter == _map.end())
  {
    return dummy;
  }
  return iter->second;
}

std::vector<RawCluster*> SvtxTrackCaloClusterMap_v1::insert(SvtxTrack* track, std::vector<RawCluster*> clusters)
{
  auto iter = _map.find(track);
  if (iter != _map.end())
  {
    iter->second = clusters;
  }
  else
  {
    _map.insert(std::make_pair(track, clusters));
  }

  return _map[track];
}

std::vector<RawCluster*> SvtxTrackCaloClusterMap_v1::insert(SvtxTrack* track, RawCluster* clus)
{
  auto iter = _map.find(track);
  if (iter != _map.end())
  {
    (iter->second).push_back(clus);
  }
  else
  {
    std::vector<RawCluster*> dummy;
    dummy.push_back(clus);
    _map.insert(std::make_pair(track, dummy));
  }

  return _map[track];
}

size_t SvtxTrackCaloClusterMap_v1::erase(SvtxTrack* track)
{
  const auto iter = _map.find(track);
  if (iter != _map.end())
  {
    _map.erase(iter);
    return 1;
  }
  else
  {
    return 0;
  }
}
