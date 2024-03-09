#ifndef FUN4ALLRAW_GL1PACKET_H
#define FUN4ALLRAW_GL1PACKET_H

#include "OfflinePacketv1.h"

#include <limits>

class  Gl1Packet: public OfflinePacketv1
  {

public:
    Gl1Packet() = default;
    ~Gl1Packet() override = default;
    using OfflinePacket::FillFrom;
    virtual void setBunchNumber(const char /*bn*/) {return;}
  virtual char getBunchNumber() const {return std::numeric_limits<char>::min();}
  virtual void FillFrom(const Gl1Packet */*pkt*/) {return;}

private:
  ClassDefOverride(Gl1Packet,1)
};

#endif
