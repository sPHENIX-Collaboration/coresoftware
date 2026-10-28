#include "MinimumBiasInfov2.h"

void MinimumBiasInfov2::identify(std::ostream& os) const
{
  os << "MinimumBiasInfo: " << std::endl;
  os << "  IsMinBias = " << (_isMinBias ? "Yes" : "No") << std::endl;
  os << "  BunchClass = " << _bunchClass << std::endl;
  os << "  T0Corrected = " << _t0Corrected << std::endl;

  return;
}

void MinimumBiasInfov2::Reset()
{
  _isMinBias = false;
  _bunchClass = MAIN;
  _t0Corrected = std::numeric_limits<float>::quiet_NaN();
}

void MinimumBiasInfov2::CopyTo(MinimumBiasInfo* mbinfo)
{
  mbinfo->setIsAuAuMinimumBias(isAuAuMinimumBias());
  mbinfo->setBunchClass(getBunchClass());
  mbinfo->setT0Corrected(getT0Corrected());
}
