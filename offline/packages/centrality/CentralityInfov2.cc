#include "CentralityInfov2.h"

CentralityInfov2::CentralityInfov2()
{
  _isMinBias = false;
}

void CentralityInfov2::identify(std::ostream& os) const
{
  os << "CentralityInfo: " << std::endl;
  os << "IsMinBias: " << (_isMinBias ? "yes":"no") << std::endl;
  return;
}

