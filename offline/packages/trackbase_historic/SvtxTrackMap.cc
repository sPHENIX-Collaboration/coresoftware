#include "SvtxTrackMap.h"

namespace 
{
  SvtxTrackMap::TrackMap DummyTrackMap;
}

SvtxTrackMap::ConstIter SvtxTrackMap::begin() const
{
  return DummyTrackMap.end();
}

SvtxTrackMap::ConstIter SvtxTrackMap::find(unsigned int /*idkey*/) const
{
  return DummyTrackMap.end();
}

SvtxTrackMap::ConstIter SvtxTrackMap::end() const
{
  return DummyTrackMap.end();
}


SvtxTrackMap::Iter SvtxTrackMap::begin()
{
  return DummyTrackMap.end();
}

SvtxTrackMap::Iter SvtxTrackMap::find(unsigned int /*idkey*/)
{
  return DummyTrackMap.end();
}

SvtxTrackMap::Iter SvtxTrackMap::end()
{
  return DummyTrackMap.end();
}
