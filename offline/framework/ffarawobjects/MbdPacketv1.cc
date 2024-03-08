#include "MbdPacketv1.h"

void MbdPacketv1::Reset()
{
  OfflinePacketv1::Reset();
  PacketEvtSequence = 0;
  femclock.fill(std::numeric_limits<uint32_t>::max());
  for (auto &row :  samples)
  {
    row.fill(std::numeric_limits<uint32_t>::max());
  }  
  return;
}

void MbdPacketv1::identify(std::ostream &os) const
{
  os << "MbdPacketv1: " << std::endl;
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
