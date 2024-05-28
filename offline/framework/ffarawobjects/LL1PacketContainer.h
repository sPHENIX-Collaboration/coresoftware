#ifndef FUN4ALLRAW_LL1PACKETCONTAINER_H
#define FUN4ALLRAW_LL1PACKETCONTAINER_H

#include <phool/PHObject.h>

class LL1Packet;

class LL1PacketContainer : public PHObject
{
 public:
  LL1PacketContainer() = default;
  virtual ~LL1PacketContainer() = default;

  virtual LL1Packet *AddPacket() { return nullptr; }
  virtual LL1Packet *AddPacket(LL1Packet *) { return nullptr; }
  virtual unsigned int get_npackets() { return 0; }
  virtual LL1Packet *getPacket(unsigned int) { return nullptr; }
  virtual LL1Packet *getPacketbyId(int) { return nullptr; }

 private:
  ClassDefOverride(LL1PacketContainer, 1)
};

#endif
