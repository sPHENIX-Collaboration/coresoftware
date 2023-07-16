#ifndef FUN4ALLRAW_PACKETLIST_H
#define FUN4ALLRAW_PACKETLIST_H

#include <iostream>
#include <map>
#include <vector>

class Packet;

class PacketList
{
 public:
  PacketList() = default;
  virtual ~PacketList();
  std::vector<Packet *>::const_iterator begin(const int packet_id);
  std::vector<Packet *>::const_iterator end(const int packet_id);
  void AddPacket(const int packet_id, Packet *pkt);
  void Reset();
  void identify(std::ostream &os = std::cout);

 private:
  std::map<int, std::vector<Packet *>> m_PacketMap;
};

#endif
