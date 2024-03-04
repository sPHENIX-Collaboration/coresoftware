#ifndef FUN4ALLRAW_MBDPACKETV1_H
#define FUN4ALLRAW_MBDPACKETV1_H

#include "MbdPacket.h"

#include <array>
#include <limits>

class  MbdPacketv1: public MbdPacket
  {

public:
    MbdPacketv1() = default;
    MbdPacketv1(MbdPacket *mbd);
    ~MbdPacketv1() override = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;

  void setFemClock(int i, uint32_t clk) override {femclock.at(i) = clk;}
  uint32_t getFemClock(int i) const override {return femclock.at(i);}
  void setSample(int ipmt, int isamp, uint32_t val) override {samples.at(isamp).at(ipmt) = val;}
  uint32_t getSample(int ipmt, int isamp) const override {return samples.at(isamp).at(ipmt);}
  void setPacketEvtSequence(int i) override {PacketEvtSequence = i;}
  int getPacketEvtSequence() const override {return PacketEvtSequence;}



  protected:
  int PacketEvtSequence {0};
  std::array<uint32_t,2> femclock {std::numeric_limits<uint32_t>::max()};
  std::array<std::array<uint32_t,128>, 31> samples{{{std::numeric_limits<uint32_t>::max()}}};

private:
  ClassDefOverride(MbdPacketv1,1)
};

#endif
