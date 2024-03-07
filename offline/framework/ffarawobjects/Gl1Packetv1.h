#ifndef FUN4ALLRAW_GL1PACKETV1_H
#define FUN4ALLRAW_GL1PACKETV1_H

#include "Gl1Packet.h"

#include <limits>

class  Gl1Packetv1: public Gl1Packet
  {

public:
    Gl1Packetv1() = default;
    ~Gl1Packetv1() override = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;
  void FillFrom(const Gl1Packet *pkt) override;

  void setBunchNumber(const char bn) override {bunchnumber = bn;}
  char getBunchNumber() const override {return bunchnumber;}

  protected:
  char bunchnumber {std::numeric_limits<char>::min()};
private:
  ClassDefOverride(Gl1Packetv1,1)
};

#endif
