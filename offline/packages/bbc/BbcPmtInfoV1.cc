#include "BbcPmtInfoV1.h"
//#include <iostream>

void BbcPmtInfoV1::Reset()
{
  Clear();
}

void BbcPmtInfoV1::Clear(Option_t* )
{
  //std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  bq = NAN;
  btt = NAN;
  btq = NAN;
}

void BbcPmtInfoV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a BbcPmtInfoV1 object" << std::endl;
  out << "Pmt: " << bpmt << ", q: " << bq << ", btt: "
      << btt << ", btq: " << btq << std::endl;
}
