#ifndef FUN4ALLRAW_PACKETLIST_H
#define FUN4ALLRAW_PACKETLIST_H

#include <iostream>
#include <set>
#include <vector>

class Packet;

class PacketList
{
 public:
  PacketList() = default;
  virtual ~PacketList();
  //  void identify(std::ostream &os = std::cout);

  // private:
  std::set<uint64_t> m_BeamClockSet;
  std::vector<Packet *> m_PacketVector;
};

#endif
