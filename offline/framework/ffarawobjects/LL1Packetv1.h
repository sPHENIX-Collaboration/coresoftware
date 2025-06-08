#ifndef FUN4ALLRAW_LL1PACKETV1_H
#define FUN4ALLRAW_LL1PACKETV1_H

#include "LL1Packet.h"

#include <array>
#include <limits>

static const int MAX_NUM_CHANNELS = 672;
static const int MAX_NUM_SAMPLES = 5;

class LL1Packetv1 : public LL1Packet
{
 public:
  LL1Packetv1();
  ~LL1Packetv1() override = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;

  int getMaxNumChannels() const override { return MAX_NUM_CHANNELS; }
  int getMaxNumSamples() const override { return MAX_NUM_SAMPLES; }

  void setNrChannels(int i) override { NrChannels = i; }
  int getNrChannels() const override { return NrChannels; }
  void setNrSamples(int i) override { NrSamples = i; }
  int getNrSamples() const override { return NrSamples; }
  void setTriggerWords(int i) override { TriggerWords = i; }
  int getTriggerWords() const override { return TriggerWords; }
  void setSlotNr(int i) override { SlotNr = i; }
  int getSlotNr() const override { return SlotNr; }
  void setCardNr(int i) override { CardNr = i; }
  int getCardNr() const override { return CardNr; }
  void setMonitor(int i) override { Monitor = i; }
  int getMonitor() const override { return Monitor; }
  void setFemWords(int i) override { FemWords = i; }
  int getFemWords() const override { return FemWords; }
  void setSums(int i) override { Sums = i; }
  int getSums() const override { return Sums; }
  void setFibers(int i) override { Fibers = i; }
  int getFibers() const override { return Fibers; }

  void setSample(int ipmt, int isamp, uint32_t val) override { samples.at(isamp).at(ipmt) = val; }
  uint32_t getSample(int ipmt, int isamp) const override { return samples.at(isamp).at(ipmt); }
  void setPacketEvtSequence(int i) override { PacketEvtSequence = i; }
  int getPacketEvtSequence() const override { return PacketEvtSequence; }
  int iValue(const int i, const std::string &what) const override;
  int iValue(const int channel, const int sample) const override;
  void dump(std::ostream &os = std::cout) const override;
  void dump_idll1_mbd(std::ostream &os = std::cout) const;
  void dump_idll1_emcal_mon3(std::ostream &os = std::cout) const;
  void dump_idll1_jet_emcal_mon1(std::ostream &os = std::cout) const;

 protected:
  int PacketEvtSequence{0};
  int NrChannels{0};
  int NrSamples{0};
  int TriggerWords{0};
  int SlotNr{0};
  int CardNr{0};
  int Monitor{0};
  int FemWords{0};
  int Sums{0};
  int Fibers{0};

  std::array<std::array<uint32_t, MAX_NUM_CHANNELS>, MAX_NUM_SAMPLES> samples{};

 private:
  ClassDefOverride(LL1Packetv1, 2)
};

#endif
