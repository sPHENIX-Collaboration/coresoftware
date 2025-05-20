#ifndef FUN4ALLRAW_LL1PACKET_H
#define FUN4ALLRAW_LL1PACKET_H

#include "OfflinePacketv1.h"

#include <limits>

class LL1Packet : public OfflinePacketv1
{
 public:
  LL1Packet() = default;
  LL1Packet(LL1Packet *pkt)
    : OfflinePacketv1(pkt)
  {
  }
  ~LL1Packet() override = default;

  virtual int getMaxNumChannels() const { return 0; }
  virtual int getMaxNumSamples() const { return 0; }
  virtual int getMaxNumModules() const { return 0; }

  virtual void setSample(int /*ipmt*/, int /*ichan*/, uint32_t /*val*/) { return; }
  virtual uint32_t getSample(int /*ipmt*/, int /*ichan*/) const { return std::numeric_limits<uint32_t>::max(); }
  virtual void setPacketEvtSequence(int /*i*/) { return; }
  virtual int getPacketEvtSequence() const { return std::numeric_limits<int>::max(); }
  virtual void setNrChannels(int /*i*/) { return; }
  virtual int getNrChannels() const { return 0; }
  virtual void setNrSamples(int /*i*/) { return; }
  virtual int getNrSamples() const { return 0; }
  virtual void setTriggerWords(int /*i*/) { return; }
  virtual int getTriggerWords() const { return 0; }
  virtual void setSlotNr(int /*i*/) { return; }
  virtual int getSlotNr() const { return 0; }
  virtual void setCardNr(int /*i*/) { return; }
  virtual int getCardNr() const { return 0; }
  virtual void setMonitor(int /*i*/) { return; }
  virtual int getMonitor() const { return 0; }
  virtual void setFemWords(int /*i*/) { return; }
  virtual int getFemWords() const { return 0; }
  virtual void setSums(int /*i*/) { return; }
  virtual int getSums() const { return 0; }
  virtual void setFibers(int /*i*/) { return; }
  virtual int getFibers() const { return 0; }

 private:
  ClassDefOverride(LL1Packet, 2)
};

#endif
