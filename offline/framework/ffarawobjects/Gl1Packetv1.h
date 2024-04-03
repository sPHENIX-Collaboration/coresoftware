#ifndef FUN4ALLRAW_GL1PACKETV1_H
#define FUN4ALLRAW_GL1PACKETV1_H

#include "Gl1Packet.h"

#include <array>
#include <limits>

class Gl1Packetv1 : public Gl1Packet
{
 public:
  Gl1Packetv1() = default;
  ~Gl1Packetv1() override = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;
  void FillFrom(const Gl1Packet *pkt) override;

  void setPacketNumber(const unsigned int i) override { packet_nr = i; }
  unsigned int getPacketNumber() const override { return packet_nr; }

  void setBunchNumber(const uint64_t i) override { BunchNumber = i; }
  uint64_t getBunchNumber() const override { return BunchNumber; }
  void setTriggerInput(const uint64_t i) override { TriggerInput = i; }
  uint64_t getTriggerInput() const override { return TriggerInput; }
  void setTriggerVector(const uint64_t i) override { TriggerVector = i; }
  uint64_t getTriggerVector() const override { return TriggerVector; }
  void setScaler(int iscal, int index, uint64_t lval) override { scaler.at(iscal).at(index) = lval; }

  int iValue(const int i) const override;
  long long lValue(const int /*i*/, const std::string &what) const override;
  long long lValue(const int i, const int j) const override;

  void dump(std::ostream &os = std::cout) const override;

 protected:
  unsigned int packet_nr{0};
  uint64_t BunchNumber{std::numeric_limits<uint64_t>::max()};
  uint64_t TriggerInput{0};
  uint64_t TriggerVector{0};
  std::array<std::array<uint64_t, 3>, 64> scaler{{{0}}};

 private:
  ClassDefOverride(Gl1Packetv1, 1)
};

#endif
