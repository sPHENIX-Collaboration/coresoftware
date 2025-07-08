#include "OfflinePacketv1.h"

OfflinePacketv1::OfflinePacketv1(OfflinePacket *pkt)
{
  setEvtSequence(pkt->getEvtSequence());
  setIdentifier(pkt->getIdentifier());
  setHitFormat(pkt->getHitFormat());
  setBCO(pkt->getBCO());
  setStatus(pkt->getStatus());
}

void OfflinePacketv1::Reset()
{
  evtseq = std::numeric_limits<int>::min();
  packetid = std::numeric_limits<int>::min();
  hitformat = std::numeric_limits<int>::min();
  bco = std::numeric_limits<uint64_t>::max();
  status = 0;
}

void OfflinePacketv1::identify(std::ostream &os) const
{
  os << "Id: " << getIdentifier() << std::endl;
  os << "Status: " << getStatus();
  if (getStatus())
  {
    std::cout << " --> Bad";
  }
  else
  {
    std::cout << " --> Good";
  }
  std::cout << std::endl;
  os << "EvtSeq: " << getEvtSequence() << std::endl;
  os << "BCO: 0x" << std::hex << getBCO() << std::dec << std::endl;
  return;
}

void OfflinePacketv1::FillFrom(const OfflinePacket *pkt)
{
  setIdentifier(pkt->getIdentifier());
  setEvtSequence(pkt->getEvtSequence());
  setHitFormat(pkt->getHitFormat());
  setBCO(pkt->getBCO());
  setStatus(pkt->getStatus());
}
