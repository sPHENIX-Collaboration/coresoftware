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

  PacketMap() = default;
  virtual ~PacketMap();
  ConstRange begin_end(const int packet_id);
  void AddPacket(const int packet_id, Packet *pkt);
  void Reset();
  void identify(std::ostream &os = std::cout);

 private:
  std::map<int, PacketList> m_PacketMap;
};

#endif
