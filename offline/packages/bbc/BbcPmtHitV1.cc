#include "BbcPmtHitV1.h"

#include <iostream>

BbcPmtHitV1::BbcPmtHitV1(const short ipmt, const float iadc, const float itdc0, const float itdc1)
  : pmt(ipmt)
  , adc(iadc)
  , tdc0(itdc0)
  , tdc1(itdc1)
{
}

void BbcPmtHitV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a BbcPmtHitV1 object" << std::endl;
  out << "Pmt: " << pmt << ", adc: " << adc << ", tdc0: "
      << tdc0 << ", tdc1: " << tdc1 << std::endl;
}
