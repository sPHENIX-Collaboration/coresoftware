#include "OfflinePacketv1.h"

OfflinePacketv1::OfflinePacketv1(OfflinePacket *pkt)
{
  setEvtSequence(pkt->getEvtSequence());
  setIdentifier(pkt->getIdentifier());
  setBCO(pkt->getBCO());
}



void OfflinePacketv1::Reset()
{
  evtseq = std::numeric_limits<int>::min();
  packetid = std::numeric_limits<int>::min();
  bco = std::numeric_limits<uint64_t>::max();
}
