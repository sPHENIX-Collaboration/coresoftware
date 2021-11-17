#include "SvtxVertex.h"

SvtxVertex::TrackSet DummyTrackSet;

SvtxVertex::ConstTrackIter SvtxVertex::begin_tracks() const
{
  return DummyTrackSet.end();
}

SvtxVertex::ConstTrackIter SvtxVertex::find_track(unsigned int /*trackid*/) const
{
  return DummyTrackSet.end();
}

SvtxVertex::ConstTrackIter SvtxVertex::end_tracks() const
{
  return DummyTrackSet.end();
}

SvtxVertex::TrackIter SvtxVertex::begin_tracks()
{
  return DummyTrackSet.end();
}

SvtxVertex::TrackIter SvtxVertex::find_track(unsigned int /*trackid*/)
{
  return DummyTrackSet.end();
}

SvtxVertex::TrackIter SvtxVertex::end_tracks()
{
  return DummyTrackSet.end();
}
