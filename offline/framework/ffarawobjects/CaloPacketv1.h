#ifndef FUN4ALLRAW_CALOPACKETV1_H
#define FUN4ALLRAW_CALOPACKETV1_H

#include "CaloPacket.h"

#include <array>
#include <limits>

static const int MAX_NUM_CHANNELS = 256;
static const int MAX_NUM_MODULES = 4;
static const int MAX_NUM_SAMPLES = 31;

class CaloPacketv1 : public CaloPacket
{
 public:
  CaloPacketv1();
  ~CaloPacketv1() override = default;

  void Reset() override;
  void identify(std::ostream &os = std::cout) const override;

  int getMaxNumChannels() const override {return MAX_NUM_CHANNELS;}
  int getMaxNumSamples() const override {return MAX_NUM_SAMPLES;}
  int getMaxNumModules() const override {return MAX_NUM_MODULES;}

  void setFemClock(int i, uint32_t clk) override { femclock.at(i) = clk; }
  uint32_t getFemClock(int i) const override { return femclock.at(i); }
  void setFemEvtSequence(int i, int evtno) override { femevt.at(i) = evtno; }
  int getFemEvtSequence(int i) const override { return femevt.at(i); }
  void setFemSlot(int i, int islot) override { femslot.at(i) = islot; }
  int getFemSlot(int i) const override { return femslot.at(i); }
  void setChecksumLsb(int i, int ival) override {checksumlsb.at(i) = ival; }
  int getChecksumLsb(int i) const override {return checksumlsb.at(i); }
  void setChecksumMsb(int i, int ival) override {checksummsb.at(i) = ival; }
  int getChecksumMsb(int i) const override {return checksummsb.at(i); }

  void setNrChannels(int i) override { NrChannels = i; }
  int getNrChannels() const override { return NrChannels; }
  void setNrSamples(int i) override { NrSamples = i; }
  int getNrSamples() const override { return NrSamples; }
  void setNrModules(int i) override { NrModules = i; }
  int getNrModules() const override { return NrModules; }
  void setEvenChecksum(int i) override { event_checksum = i; }
  int getEvenChecksum() const override { return event_checksum; }
  void setOddChecksum(int i) override { odd_checksum = i; }
  int getOddChecksum() const override { return odd_checksum; }
  void setCalcEvenChecksum(int i) override { calc_event_checksum = i; }
  int getCalcEvenChecksum() const override { return calc_event_checksum; }
  void setCalcOddChecksum(int i) override { calc_odd_checksum = i; }
  int getCalcOddChecksum() const override { return calc_odd_checksum; }
  void setModuleAddress(int i) override { module_address = i; }
  int getModuleAddress() const override { return module_address; }
  void setDetId(int i) override { detid = i; }
  int getDetId() const override { return detid; }
  bool getSuppressed(int channel) const override {return isZeroSuppressed.at(channel);}
  void setSuppressed(int channel, bool bb) override {isZeroSuppressed.at(channel) = bb;}
  void setPre(int channel, uint32_t ival) override {pre.at(channel) = ival;}
  uint32_t getPre(int channel) const override {return pre.at(channel);}
  void setPost(int channel, uint32_t ival) override {post.at(channel) = ival;}
  uint32_t getPost(int channel) const override {return post.at(channel);}

  void setSample(int ipmt, int isamp, uint32_t val) override { samples.at(isamp).at(ipmt) = val; }
  uint32_t getSample(int ipmt, int isamp) const override { return samples.at(isamp).at(ipmt); }
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
  int NrModules{0};
  int event_checksum{0};
  int odd_checksum{0};
  int calc_event_checksum{0};
  int calc_odd_checksum{0};
  int module_address{0};
  int detid{0};

  std::array<uint32_t, MAX_NUM_MODULES> femclock{};
  std::array<uint32_t, MAX_NUM_MODULES> femevt {};
  std::array<uint32_t, MAX_NUM_MODULES> femslot {};
  std::array<uint32_t, MAX_NUM_MODULES> checksumlsb {};
  std::array<uint32_t, MAX_NUM_MODULES> checksummsb {};

  std::array<std::array<uint32_t, MAX_NUM_CHANNELS>, MAX_NUM_SAMPLES> samples {};
  std::array<bool,MAX_NUM_CHANNELS> isZeroSuppressed {};
  std::array<uint32_t, MAX_NUM_CHANNELS> pre {};
  std::array<uint32_t, MAX_NUM_CHANNELS> post {};

 private:
  ClassDefOverride(CaloPacketv1, 1)
};

#endif
