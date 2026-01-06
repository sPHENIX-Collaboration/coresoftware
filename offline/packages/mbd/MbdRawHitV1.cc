#include "MbdRawHitV1.h"

void MbdRawHitV1::Reset()
{
  Clear();
}

void MbdRawHitV1::Clear(Option_t* /*unused*/)
{
  std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  badc = std::numeric_limits<float>::quiet_NaN();
  bttdc = std::numeric_limits<float>::quiet_NaN();
  bqtdc = std::numeric_limits<float>::quiet_NaN();
}

void MbdRawHitV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a MbdRawHitV1 object" << std::endl;
  out << "Pmt: " << bpmt << ", adc: " << badc << ", ttdc: "
      << bttdc << ", bqtdc: " << bqtdc << std::endl;
}
