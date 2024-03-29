#include "SvtxTrackCaloClusterMap.h"

namespace
{
  SvtxTrackCaloClusterMap::Map DummyTrackMap;
}

SvtxTrackCaloClusterMap::ConstIter SvtxTrackCaloClusterMap::begin() const
{
  return DummyTrackMap.end();
}

SvtxTrackCaloClusterMap::ConstIter SvtxTrackCaloClusterMap::find(SvtxTrack* /*unused*/) const
{
  return DummyTrackMap.end();
}

SvtxTrackCaloClusterMap::ConstIter SvtxTrackCaloClusterMap::end() const
{
  return DummyTrackMap.end();
}

SvtxTrackCaloClusterMap::Iter SvtxTrackCaloClusterMap::begin()
{
  return DummyTrackMap.end();
}

SvtxTrackCaloClusterMap::Iter SvtxTrackCaloClusterMap::find(SvtxTrack* /*unused*/)
{
  return DummyTrackMap.end();
}

SvtxTrackCaloClusterMap::Iter SvtxTrackCaloClusterMap::end()
{
  return DummyTrackMap.end();
}
