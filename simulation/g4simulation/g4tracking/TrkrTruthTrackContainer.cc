#include "TrkrTruthTrackContainer.h"

namespace
{
  TrkrTruthTrackContainer::Map dummy_map;
}

TrkrTruthTrackContainer::ConstRange TrkrTruthTrackContainer::getTruthTrackRange() const
{
  return std::make_pair(dummy_map.begin(), dummy_map.end());
}

TrkrTruthTrackContainer::Map& TrkrTruthTrackContainer::getMap()
{
  return dummy_map;
}

int TrkrTruthTrackContainer::nhw_cc() { return 106; };
