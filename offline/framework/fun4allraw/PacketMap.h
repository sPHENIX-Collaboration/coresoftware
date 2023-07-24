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
  typedef std::vector<Packet *>::const_iterator ConstIterator;
  typedef std::pair<ConstIterator, ConstIterator> ConstRange;
  typedef std::map<int, PacketList>::const_iterator PacketIterator;
  typedef std::pair<PacketIterator, PacketIterator> PacketRange;

  PacketMap() = default;
  virtual ~PacketMap();
  ConstRange begin_end(const int packet_id);
  PacketRange first_last_packet();
  void AddPacket(const int packet_id, Packet *pkt);
  void AddBclk(const int packet_id, uint64_t bclk);
  //  uint64_t GetCurrentBeamClock();
  void Reset();
  void identify(std::ostream &os = std::cout);

 private:
  std::map<int, PacketList> m_PacketMap;
};

#endif
