#include "TrkrTruthTrack.h"

namespace
{
  std::vector<TrkrDefs::cluskey> dummyClusters;
}

std::vector<TrkrDefs::cluskey>& TrkrTruthTrack::getClusters()
{
  return dummyClusters;
}
