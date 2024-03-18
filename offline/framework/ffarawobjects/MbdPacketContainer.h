#ifndef FUN4ALLRAW_MBDRAWHITCONTAINER_H
#define FUN4ALLRAW_MBDRAWHITCONTAINER_H

#include <phool/PHObject.h>

class MbdPacket;

class MbdPacketContainer : public PHObject
{
 public:
  MbdPacketContainer() = default;
  virtual ~MbdPacketContainer() = default;

  virtual MbdPacket *AddPacket() { return nullptr; }
  virtual MbdPacket *AddPacket(MbdPacket *) { return nullptr; }
  virtual unsigned int get_npackets() { return 0; }
  virtual MbdPacket *getPacket(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(MbdPacketContainer, 1)
};

#endif
