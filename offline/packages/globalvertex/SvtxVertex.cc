#include "SvtxVertex.h"

SvtxVertex::TrackSet trackSet;

SvtxVertex::ConstTrackIter SvtxVertex::begin_tracks() const
{
  return trackSet.end();
}

SvtxVertex::ConstTrackIter SvtxVertex::find_track(unsigned int /*trackid*/) const
{
  return trackSet.end();
}

SvtxVertex::ConstTrackIter SvtxVertex::end_tracks() const
{
  return trackSet.end();
}

SvtxVertex::TrackIter SvtxVertex::begin_tracks()
{
  return trackSet.end();
}

SvtxVertex::TrackIter SvtxVertex::find_track(unsigned int /*trackid*/)
{
  return trackSet.end();
}

SvtxVertex::TrackIter SvtxVertex::end_tracks()
{
  return trackSet.end();
}
