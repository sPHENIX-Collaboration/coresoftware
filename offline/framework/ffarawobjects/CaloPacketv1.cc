#include "CaloPacketv1.h"

void CaloPacketv1::Reset()
{
  OfflinePacketv1::Reset();
  femclock.fill(std::numeric_limits<uint32_t>::max());
  return;
}

void CaloPacketv1::identify(std::ostream &os) const
{
  os << "CaloPacketv1: " << std::endl;
    OfflinePacketv1::identify(os);
  os << "Pkt Event no: " <<  getPacketEvtSequence() << std::endl;
  os << "FEM clk: " << std::hex << getFemClock(0) << ", "
     << getFemClock(1) << std::dec << std::endl;
/*
  for (auto &iter :  samples)
  {
    for (auto &iter2 : iter)
    {
    std::cout << "sample: " << iter2 << std::endl;
    }
  }
*/
}
