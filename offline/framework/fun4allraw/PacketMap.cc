#include "PacketMap.h"

#include <Event/packet.h>

std::vector<Packet *> empty_vector;
std::set<uint64_t> empty_set;

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

void PacketMap::AddBclk(const int packet_id, uint64_t bclk)
{
  m_PacketMap[packet_id].m_BeamClockSet.insert(bclk);
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
    for (auto const setiter : iditer.second.m_BeamClockSet)
    {
      os << "BeamClock 0x" << std::hex << setiter << std::dec << std::endl;
    }
  }
}

PacketMap::PacketRange PacketMap::begin_end(const int packet_id)
{
  auto mapiter = m_PacketMap.find(packet_id);
  if (mapiter == m_PacketMap.end())
  {
    return std::make_pair(empty_vector.end(), empty_vector.end());
  }
  return std::make_pair(mapiter->second.m_PacketVector.begin(), mapiter->second.m_PacketVector.end());
}

PacketMap::BclkRange PacketMap::begin_end_bclk(const int packet_id)
{
  auto mapiter = m_PacketMap.find(packet_id);
  if (mapiter == m_PacketMap.end())
  {
    return std::make_pair(empty_set.end(), empty_set.end());
  }
  return std::make_pair(mapiter->second.m_BeamClockSet.begin(), mapiter->second.m_BeamClockSet.end());
}

PacketMap::PacketListRange PacketMap::first_last_packet()
{
  return std::make_pair(m_PacketMap.begin(), m_PacketMap.end());
}
