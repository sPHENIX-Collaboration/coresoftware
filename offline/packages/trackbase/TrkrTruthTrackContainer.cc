#include "TrkrTruthTrackContainer.h"

namespace
{
  TrkrTruthTrackContainer::Vector dummy_vector;
}

TrkrTruthTrackContainer::ConstRange TrkrTruthTrackContainer::getTruthTrackRange() const
{
  return std::make_pair(dummy_vector.begin(), dummy_vector.end());
}

TrkrTruthTrackContainer::Vector& TrkrTruthTrackContainer::getTruthTracks()
{
  return dummy_vector;
}
