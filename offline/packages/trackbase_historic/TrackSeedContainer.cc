#include "TrackSeedContainer.h"

TrackSeedContainer::Container container;

TrackSeedContainer::ConstIter TrackSeedContainer::begin() const
{
  return container.end();
}

TrackSeedContainer::ConstIter TrackSeedContainer::find(const std::size_t /*unused*/) const
{
  return container.end();
}

TrackSeedContainer::ConstIter TrackSeedContainer::end() const
{
  return container.end();
}

TrackSeedContainer::Iter TrackSeedContainer::begin()
{
  return container.end();
}

TrackSeedContainer::Iter TrackSeedContainer::find(const std::size_t /*unused*/)
{
  return container.end();
}

TrackSeedContainer::Iter TrackSeedContainer::end()
{
  return container.end();
}
