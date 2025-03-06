#include "BbcPmtInfoV1.h"

void BbcPmtInfoV1::Reset()
{
  Clear();
}

void BbcPmtInfoV1::Clear(Option_t* /*unused*/)
{
  // std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  bq = std::numeric_limits<Float_t>::quiet_NaN();
  btt = std::numeric_limits<Float_t>::quiet_NaN();
  btq = std::numeric_limits<Float_t>::quiet_NaN();
}

void BbcPmtInfoV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a BbcPmtInfoV1 object" << std::endl;
  out << "Pmt: " << bpmt << ", q: " << bq << ", btt: "
      << btt << ", btq: " << btq << std::endl;
}
