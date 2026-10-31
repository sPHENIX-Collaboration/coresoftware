#ifndef PHGARFIELDMAPFINGERPRINT_H
#define PHGARFIELDMAPFINGERPRINT_H

#include "PHGarfieldAdaptiveMap.h"

#include <string>
#include <string_view>

namespace PHGarfieldMapFingerprint
{
  // Identical primitive order for native and capacity-independent identities.
  std::string Adaptive(const PHGarfieldAdaptiveMap::Data &, std::string_view policy_prefix);
}  // namespace PHGarfieldMapFingerprint

#endif
