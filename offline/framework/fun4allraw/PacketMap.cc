#include "PacketMap.h"

#include <Event/packet.h>

std::vector<Packet *> empty_vector;

PacketMap::~PacketMap()
{
  // for (auto iditer : m_PacketMap)
  // {
  //   for (auto pktiter : iditer.second)
  //   {
  //     delete pktiter;
  //   }
  //   iditer.second.clear();
  // }
  m_PacketMap.clear();
}

void PacketMap::AddPacket(const int packet_id, Packet *pkt)
{
  m_PacketMap[packet_id].m_PacketVector.push_back(pkt);
  return;
}

void PacketMap::Reset()
{
//   for (auto iditer : m_PacketMap)
//   {
// //    for (auto pktiter : iditer.second)
//     {
// //      delete pktiter;
//     }
//     iditer.second.clear();
//   }
  m_PacketMap.clear();
}

void PacketMap::identify(std::ostream &os)
{
  for (auto const &iditer : m_PacketMap)
  {
    os << "Packet " << iditer.first << " number of packets: "
       << iditer.second.m_PacketVector.size() << std::endl;
  }
}

PacketMap::ConstRange PacketMap::begin_end(const int packet_id)
{
  auto mapiter = m_PacketMap.find(packet_id);
  if (mapiter == m_PacketMap.end())
  {
    return std::make_pair(empty_vector.end(), empty_vector.end());
  }
  return std::make_pair(mapiter->second.m_PacketVector.begin(),mapiter->second.m_PacketVector.end());
}
