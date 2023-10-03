#include "CentralityInfov2.h"

void CentralityInfov2::identify(std::ostream& os) const
{
  os << "  CentralityInfo: " << std::endl;
  os << "      IsMinBias: " << (_isMinBias ? "yes" : "no") << std::endl;
  os << "      Centrile: " << (has_centile(CentralityInfo::PROP::mbd_NS) ? get_centile(CentralityInfo::PROP::mbd_NS) : -999.99) << std::endl;

  return;
}
