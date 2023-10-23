#include "MbdPmtHitV1.h"


void MbdPmtHitV1::Reset()
{
  Clear();
}

void MbdPmtHitV1::Clear(Option_t* )
{
  //std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  bq = std::numeric_limits<float>::quiet_NaN();
  btt = std::numeric_limits<float>::quiet_NaN();
  btq = std::numeric_limits<float>::quiet_NaN();
}


void MbdPmtHitV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a MbdPmtHitV1 object" << std::endl;
  out << "Pmt: " << bpmt << ", Q: " << bq << ", tt: "
      << btt << ", btq: " << btq << std::endl;
}
