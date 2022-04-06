#include "SvtxTrackSeedContainer.h"

SvtxTrackSeedContainer::TrackSeedContainer container;

SvtxTrackSeedContainer::Iter SvtxTrackSeedContainer::erase(const std::size_t)
{
  return container.end();
}

SvtxTrackSeedContainer::ConstIter SvtxTrackSeedContainer::begin() const
{
  return container.end();
}

SvtxTrackSeedContainer::ConstIter SvtxTrackSeedContainer::find(const std::size_t) const
{
  return container.end();
}

SvtxTrackSeedContainer::ConstIter SvtxTrackSeedContainer::end() const
{
  return container.end();
}


SvtxTrackSeedContainer::Iter SvtxTrackSeedContainer::begin()
{
  return container.end();
}

SvtxTrackSeedContainer::Iter SvtxTrackSeedContainer::find(const std::size_t)
{
  return container.end();
}

SvtxTrackSeedContainer::Iter SvtxTrackSeedContainer::end()
{
  return container.end();
}
