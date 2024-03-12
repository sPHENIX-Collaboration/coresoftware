#ifndef FUN4ALLRAW_OFFLINEPACKET_H
#define FUN4ALLRAW_OFFLINEPACKET_H

#include <phool/PHObject.h>

#include <limits>

class  OfflinePacket: public PHObject
  {


public:
    OfflinePacket() = default;
  virtual ~OfflinePacket() = default;
  virtual void FillFrom(const OfflinePacket */*pkt*/) {return;}

  virtual int getIdentifier() const {return std::numeric_limits<int>::min();}
  virtual void setIdentifier(const int) {return;}
  virtual int getEvtSequence() const  {return std::numeric_limits<int>::min();}
  virtual void setEvtSequence(const int)  {return;}
  virtual uint64_t getBCO ()  const  {return std::numeric_limits<uint64_t>::max();}
  virtual void setBCO (const uint64_t) {return;}

private:
  ClassDefOverride(OfflinePacket,1)
};

#endif
