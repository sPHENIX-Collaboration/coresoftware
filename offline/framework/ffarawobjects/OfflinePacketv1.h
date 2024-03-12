#ifndef FUN4ALLRAW_OFFLINEPACKETV1_H
#define FUN4ALLRAW_OFFLINEPACKETV1_H

#include "OfflinePacket.h"

#include <limits>

class  OfflinePacketv1: public OfflinePacket
  {

public:
    OfflinePacketv1() = default;
    OfflinePacketv1(OfflinePacket *);
  ~OfflinePacketv1() override = default;
  void Reset() override;

  int getIdentifier() const override {return packetid;}
  void setIdentifier(const int i) override {packetid = i;}
  int getEvtSequence() const override {return evtseq;}
  void setEvtSequence(const int i) override {evtseq = i;}
  uint64_t getBCO()  const override {return bco;}
  void setBCO(const uint64_t ui)  override {bco = ui;}

  protected:
  int evtseq {std::numeric_limits<int>::min()};
  int packetid {std::numeric_limits<int>::min()};
  uint64_t bco {std::numeric_limits<uint64_t>::max()};

private:
  ClassDefOverride(OfflinePacketv1,1)
};

#endif
