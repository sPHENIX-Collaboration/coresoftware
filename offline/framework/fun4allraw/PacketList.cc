#include "PacketList.h"

#include <Event/packet.h>

PacketList::~PacketList()
{
  for (auto iditer: m_PacketMap)
  {
    for (auto pktiter: iditer.second)
    {
      delete pktiter;
    }
    iditer.second.clear();
  }
  m_PacketMap.clear();
}

void PacketList::AddPacket(const int packet_id, Packet *pkt)
{
  m_PacketMap[packet_id].push_back(pkt);
  return;
}

void PacketList::Reset()
{
  for (auto iditer: m_PacketMap)
  {
    for (auto pktiter: iditer.second)
    {
      delete pktiter;
    }
    iditer.second.clear();
  }
  m_PacketMap.clear();
}

void PacketList::identify(std::ostream &os)
{
  for (auto iditer: m_PacketMap)
  {
    os << "Packet " << iditer.first << " number of packets: " 
       << iditer.second.size() << std::endl;
  }
}
