#include "PacketList.h"

#include <Event/packet.h>

// std::vector<Packet *> empty_vector;

PacketList::~PacketList()
{
  // for (auto iditer : m_PacketList)
  // {
  //   for (auto pktiter : m_PacketVector)
  //   {
  //     delete pktiter;
  //   }
  //   iditer.second.clear();
  // }
  // m_PacketList.clear();
}

// void PacketList::identify(std::ostream &os)
// {
//   for (auto const &iditer : m_PacketList)
//   {
//     os << "Packet " << iditer.first << " number of packets: "
//        << iditer.second.size() << std::endl;
//   }
// }
