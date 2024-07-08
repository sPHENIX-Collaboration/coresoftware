#ifndef FUN4ALLRAW_CALOPACKETCONTAINER_H
#define FUN4ALLRAW_CALOPACKETCONTAINER_H

#include <phool/PHObject.h>

#include <limits>

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
  virtual void setEvtSequence(const int) {return;}
  virtual int getEvtSequence() const {return std::numeric_limits<int>::min();}
  virtual void setStatus(const unsigned int) {return;}
  virtual unsigned int getStatus() const {return 0;}
  virtual void deletePacketAt(int) {return;}
  virtual void deletePacket(CaloPacket *) {return;}

 private:
  ClassDefOverride(CaloPacketContainer, 1)
};

#endif
