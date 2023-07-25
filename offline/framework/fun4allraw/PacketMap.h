#ifndef FUN4ALLRAW_PACKETMAP_H
#define FUN4ALLRAW_PACKETMAP_H

#include <iostream>
#include <map>
#include <vector>

#include "PacketList.h"

class Packet;

class PacketMap
{
 public:
  typedef std::vector<Packet *>::iterator PacketIterator;
  typedef std::pair<PacketIterator, PacketIterator> PacketRange;
  typedef std::set<uint64_t>::const_iterator BclkIterator;
  typedef std::pair<BclkIterator, BclkIterator> BclkRange;
  typedef std::map<int, PacketList>::const_iterator PacketListIterator;
  typedef std::pair<PacketListIterator, PacketListIterator> PacketListRange;

  PacketMap() = default;
  virtual ~PacketMap();
  PacketRange begin_end(const int packet_id);
  BclkRange begin_end_bclk(const int packet_id);
  PacketListRange first_last_packet();
  void AddPacket(const int packet_id, Packet *pkt);
  void AddBclk(const int packet_id, uint64_t bclk);
  //  uint64_t GetCurrentBeamClock();
  void Reset();
  void identify(std::ostream &os = std::cout);

 private:
  std::map<int, PacketList> m_PacketMap;
};

#endif
