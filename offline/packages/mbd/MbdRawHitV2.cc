#include "MbdRawHitV2.h"

void MbdRawHitV2::Reset()
{
  Clear();
}

void MbdRawHitV2::Clear(Option_t* /*unused*/)
{
  //std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  fitstat = 0;
  badc = std::numeric_limits<float>::quiet_NaN();
  bttdc = std::numeric_limits<float>::quiet_NaN();
  bqtdc = std::numeric_limits<float>::quiet_NaN();
}

void MbdRawHitV2::identify(std::ostream& out) const
{
  out << "identify yourself: I am a MbdRawHitV2 object" << std::endl;
  out << "Pmt: " << bpmt << ", adc: " << badc << ", ttdc: "
      << bttdc << ", bqtdc: " << bqtdc << std::endl;
}
