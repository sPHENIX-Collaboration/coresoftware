#include "MinimumBiasInfov1.h"

#include <ostream>

void MinimumBiasInfov1::identify(std::ostream& os) const
{
  os << "MinimumBiasInfo: " << std::endl;
  os << "  IsMinBias = " << (_isMinBias ? "Yes" : "No") << std::endl;

  return;
}

void MinimumBiasInfov1::CopyTo(MinimumBiasInfo *mbinfo)
{
  mbinfo->setIsAuAuMinimumBias(isAuAuMinimumBias());
}
