#include "MbdPacketv1.h"

MbdPacketv1::MbdPacketv1(MbdPacket *pkt):
  MbdPacket(pkt)
{
  for (auto i = 0; i < 2; i++)
  {
    setFemClock(i,pkt->getFemClock(i));
  }
}


void MbdPacketv1::Reset()
{
  OfflinePacketv1::Reset();
  femclock.fill(std::numeric_limits<uint32_t>::max());
  return;
}

void MbdPacketv1::identify(std::ostream &os) const
{
  os << "MbdPacketv1: " << getIdentifier() << ", Evt sequence: " << getEvtSequence()
     << ", BCO: 0x" << std::hex <<  getBCO () << std::dec << std::endl;
  os << "Pkt Event no: " <<  getPacketEvtSequence() << std::endl;
  os << "FEM clk: " << std::hex << getFemClock(0) << ", " << getFemClock(1) << std::endl;
}
