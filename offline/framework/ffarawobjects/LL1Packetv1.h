#ifndef FUN4ALLRAW_LL1PACKETV1_H
#define FUN4ALLRAW_LL1PACKETV1_H

#include "LL1Packet.h"

#include <array>
#include <limits>

static const int MAX_NUM_CHANNELS = 256;
static const int MAX_NUM_MODULES = 4;
static const int MAX_NUM_SAMPLES = 31;

class LL1Packetv1 : public LL1Packet
{
 public:
  LL1Packetv1();
  ~LL1Packetv1() override = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;

  int getMaxNumChannels() const override {return MAX_NUM_CHANNELS;}
  int getMaxNumSamples() const override {return MAX_NUM_SAMPLES;}
  int getMaxNumModules() const override {return MAX_NUM_MODULES;}

  void setNrChannels(int i) override { NrChannels = i; }
  int getNrChannels() const override { return NrChannels; }
  void setNrSamples(int i) override { NrSamples = i; }
  int getNrSamples() const override { return NrSamples; }

//  void setSample(int ipmt, int isamp, uint32_t val) override { samples.at(isamp).at(ipmt) = val; }
//  uint32_t getSample(int ipmt, int isamp) const override { return samples.at(isamp).at(ipmt); }
  void setPacketEvtSequence(int i) override { PacketEvtSequence = i; }
  int getPacketEvtSequence() const override { return PacketEvtSequence; }
  int iValue(const int i, const std::string &what) const override;
  int iValue(const int channel, const int sample) const override;
  void dump(std::ostream &os = std::cout) const override;
  void dump93(std::ostream &os = std::cout) const;
  void dump172(std::ostream &os = std::cout) const;

 protected:
  int PacketEvtSequence{0};
  int NrChannels{0};
  int NrSamples{0};

 private:
  ClassDefOverride(LL1Packetv1, 1)
};

#endif
