#ifndef FUN4ALLRAW_MBDRAWHITCONTAINER_H
#define FUN4ALLRAW_MBDRAWHITCONTAINER_H

#include <phool/PHObject.h>

class OfflinePacket;

class  MbdPacketContainer: public PHObject
{
public:
  MbdPacketContainer() = default;
  virtual ~MbdPacketContainer() = default;

  virtual OfflinePacket *AddHit() {return nullptr;}
  virtual OfflinePacket *AddHit(OfflinePacket *) {return nullptr;}
  virtual unsigned int get_npackets() {return 0;}
  virtual OfflinePacket *get_packet(unsigned int) {return nullptr;}

private:
  ClassDefOverride(MbdPacketContainer,1)
};

#endif
