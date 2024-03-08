#ifndef FUN4ALLRAW_GL1PACKETV1_H
#define FUN4ALLRAW_GL1PACKETV1_H

#include "OfflinePacketv1.h"

#include <limits>

class  Gl1Packetv1: public OfflinePacketv1
  {

public:
    Gl1Packetv1() = default;
    ~Gl1Packetv1() override = default;

  void Reset() override;

  virtual void setBunchNumber(const char bn) {bunchnumber = bn;}
  virtual char getBunchNumber() const {return bunchnumber;}

  protected:
  char bunchnumber {std::numeric_limits<char>::min()};
private:
  ClassDefOverride(Gl1Packetv1,1)
};

#endif
