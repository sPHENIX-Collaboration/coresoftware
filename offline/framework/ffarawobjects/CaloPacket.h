#ifndef FUN4ALLRAW_CALOPACKET_H
#define FUN4ALLRAW_CALOPACKET_H

#include "OfflinePacketv1.h"

#include <array>
#include <limits>

class  CaloPacket: public OfflinePacketv1
  {

public:
    CaloPacket() = default;
  CaloPacket(CaloPacket *pkt): OfflinePacketv1(pkt) {}
    ~CaloPacket() override = default;

  virtual void setFemClock(int /*i*/, uint32_t /*clk*/) {return;}
  virtual uint32_t getFemClock(int /*i*/) const {return std::numeric_limits<uint32_t>::max();} 
  virtual void setSample(int /*ipmt*/, int /*ichan*/, uint32_t /*val*/) {return;}
  virtual uint32_t getSample(int /*ipmt*/, int /*ichan*/) const {return std::numeric_limits<uint32_t>::max();} 
  virtual void setPacketEvtSequence(int /*i*/) {return;}
  virtual int getPacketEvtSequence() const {return std::numeric_limits<int>::max();}
  virtual void setFemEvtSequence(int /*i*/,int /*i*/) {return;}
  virtual int getFemEvtSequence(int /*i*/) const {return std::numeric_limits<int>::max();}

private:
  ClassDefOverride(CaloPacket,1)
};

#endif
