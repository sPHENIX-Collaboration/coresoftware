#include "MbdPmtHitV2.h"

void MbdPmtHitV2::Reset()
{
  Clear();
}

void MbdPmtHitV2::Clear(Option_t* /*unused*/)
{
  //std::cout << "clearing " << bpmt << std::endl;
  bpmt = -1;
  fitstat = 0;
  bq = std::numeric_limits<float>::quiet_NaN();
  btt = std::numeric_limits<float>::quiet_NaN();
  btq = std::numeric_limits<float>::quiet_NaN();
}

void MbdPmtHitV2::identify(std::ostream& out) const
{
  out << "identify yourself: I am a MbdPmtHitV2 object" << std::endl;
  out << "Pmt: " << bpmt << ", fitstat 0x" << std::hex << fitstat << std::dec
      << ", Q: " << bq << ", tt: " << btt << ", btq: " << btq << std::endl;
}
