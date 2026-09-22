#ifndef TRACKRECO_PHSILICONSEEDPRUNERHELPER_H
#define TRACKRECO_PHSILICONSEEDPRUNERHELPER_H

#include <gsl/gsl_rng.h>

#include <cstddef>
#include <vector>

class TrackSeedContainer;

namespace PHSiliconSeedPrunerHelper
{
  struct Result
  {
    std::vector<std::size_t> selectedSeedIndices;
    std::vector<std::size_t> uncertifiedSeedIndices;
  };

  Result SelectSeeds(
      const TrackSeedContainer& container,
      const std::vector<std::size_t>& seedIndices,
      gsl_rng* rng,
      std::size_t mvtxLayerCount,
      std::size_t maxRepresentatives,
      std::size_t searchBudget,
      std::size_t heuristicRestarts,
      std::size_t debugMinimumGroupSize,
      bool debugVerbosity,
      bool debugDetailedVerbosity);
}  // namespace PHSiliconSeedPrunerHelper

#endif  // TRACKRECO_PHSILICONSEEDPRUNERHELPER_H
