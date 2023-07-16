#include "PacketList.h"

#include <Event/packet.h>

std::vector<Packet *> empty_vector;

PacketList::~PacketList()
{
  for (auto iditer : m_PacketMap)
  {
    for (auto pktiter : iditer.second)
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
  for (auto iditer : m_PacketMap)
  {
    for (auto pktiter : iditer.second)
    {
      delete pktiter;
    }
    iditer.second.clear();
  }
  m_PacketMap.clear();
}

void PacketList::identify(std::ostream &os)
{
  for (auto const &iditer : m_PacketMap)
  {
    os << "Packet " << iditer.first << " number of packets: "
       << iditer.second.size() << std::endl;
  }
}

std::vector<Packet *>::const_iterator PacketList::begin(const int packet_id)
{
  auto mapiter = m_PacketMap.find(packet_id);
  if (mapiter == m_PacketMap.end())
  {
    return empty_vector.begin();
  }
  return mapiter->second.begin();
}

std::vector<Packet *>::const_iterator PacketList::end(const int packet_id)
{
  auto mapiter = m_PacketMap.find(packet_id);
  if (mapiter == m_PacketMap.end())
  {
    return empty_vector.end();
  }
  return mapiter->second.end();
}
