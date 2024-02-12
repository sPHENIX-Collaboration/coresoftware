#ifndef FUN4ALLRAW_MBDPACKETV1_H
#define FUN4ALLRAW_MBDPACKETV1_H

#include "OfflinePacketv1.h"

#include <limits>

class  MbdPacketv1: public OfflinePacketv1
  {

public:
    MbdPacketv1() = default;
    ~MbdPacketv1() override = default;

  void Reset() override;

  virtual void setBunchNumber(const char bn) {bunchnumber = bn;}
  virtual char getBunchNumber() const {return bunchnumber;}

  protected:
  char bunchnumber {std::numeric_limits<char>::min()};
private:
  ClassDefOverride(MbdPacketv1,1)
};

#endif
