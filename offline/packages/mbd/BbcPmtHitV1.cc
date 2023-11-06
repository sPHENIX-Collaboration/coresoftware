#include "BbcPmtHitV1.h"


void BbcPmtHitV1::Reset()
{
  Clear();
}

void BbcPmtHitV1::Clear(Option_t* )
{
  //std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  bq = std::numeric_limits<float>::quiet_NaN();
  btt = std::numeric_limits<float>::quiet_NaN();
  btq = std::numeric_limits<float>::quiet_NaN();
}


void BbcPmtHitV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a BbcPmtHitV1 object" << std::endl;
  out << "Pmt: " << bpmt << ", Q: " << bq << ", tt: "
      << btt << ", btq: " << btq << std::endl;
}
