#ifndef FUN4ALLRAW_GL1PACKETV3_H
#define FUN4ALLRAW_GL1PACKETV3_H

#include "Gl1Packetv2.h"

#include <array>
#include <limits>

class Gl1Packetv3 : public Gl1Packetv2
{
 public:
  Gl1Packetv3() = default;
  ~Gl1Packetv3() override = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;
  void FillFrom(const Gl1Packet *pkt) override;

  void setGTMAllBusyVector(const uint64_t i) override { GTMAllBusyVector = i; }
  uint64_t getGTMAllBusyVector() const override { return GTMAllBusyVector; }

  long long lValue(const int /*i*/, const std::string &what) const override;

  void dump(std::ostream &os = std::cout) const override;

 protected:
  uint64_t GTMAllBusyVector{0};

 private:
  ClassDefOverride(Gl1Packetv3, 1)
};

#endif
