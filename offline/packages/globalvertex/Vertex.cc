#include "Vertex.h"

Vertex::TrackSet DummyTrackSet;

Vertex::ConstTrackIter Vertex::begin_tracks() const
{
  return DummyTrackSet.end();
}

Vertex::ConstTrackIter Vertex::find_track(unsigned int /*trackid*/) const
{
  return DummyTrackSet.end();
}

Vertex::ConstTrackIter Vertex::end_tracks() const
{
  return DummyTrackSet.end();
}

Vertex::TrackIter Vertex::begin_tracks()
{
  return DummyTrackSet.end();
}

Vertex::TrackIter Vertex::find_track(unsigned int /*trackid*/)
{
  return DummyTrackSet.end();
}

Vertex::TrackIter Vertex::end_tracks()
{
  return DummyTrackSet.end();
}
