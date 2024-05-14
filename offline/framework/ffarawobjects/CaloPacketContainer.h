#ifndef FUN4ALLRAW_CALORAWHITCONTAINER_H
#define FUN4ALLRAW_CALORAWHITCONTAINER_H

#include <phool/PHObject.h>

class CaloPacket;

class CaloPacketContainer : public PHObject
{
 public:
  CaloPacketContainer() = default;
  virtual ~CaloPacketContainer() = default;

  virtual CaloPacket *AddPacket() { return nullptr; }
  virtual CaloPacket *AddPacket(CaloPacket *) { return nullptr; }
  virtual unsigned int get_npackets() { return 0; }
  virtual CaloPacket *getPacket(unsigned int) { return nullptr; }
  virtual CaloPacket *getPacketbyId(int) { return nullptr; }

 private:
  ClassDefOverride(CaloPacketContainer, 1)
};

#endif
