#include "Gl1Packetv1.h"

void Gl1Packetv1::Reset()
{
  OfflinePacketv1::Reset();
  bunchnumber = std::numeric_limits<char>::min();
  return;
}

void Gl1Packetv1::identify(std::ostream &os) const
{
  os << "Gl1Packetv1: " << std::endl;
  OfflinePacketv1::identify(os);
  os << "bunch number: " << bunchnumber*1 << std::endl;
  return;
}

void Gl1Packetv1::FillFrom(const Gl1Packet *pkt)
{
  setBunchNumber(pkt->getBunchNumber());
  OfflinePacketv1::FillFrom(pkt);
}
