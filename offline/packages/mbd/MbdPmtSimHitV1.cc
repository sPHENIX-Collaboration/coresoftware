#include "MbdPmtSimHitV1.h"

void MbdPmtSimHitV1::Reset()
{
  Clear();
}

void MbdPmtSimHitV1::Clear(Option_t* /*unused*/)
{
  // std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  bq = std::numeric_limits<float>::quiet_NaN();
  btt = std::numeric_limits<float>::quiet_NaN();
  btq = std::numeric_limits<float>::quiet_NaN();
  bnpe = std::numeric_limits<float>::quiet_NaN();
}

void MbdPmtSimHitV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a MbdPmtSimHitV1 object" << std::endl;
  out << "Pmt: " << bpmt << ", Q: " << bq << ", tt: "
      << btt << ", btq: " << btq << ", bnpe: " << bnpe << std::endl;
}
