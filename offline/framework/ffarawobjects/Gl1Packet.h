#ifndef FUN4ALLRAW_GL1PACKET_H
#define FUN4ALLRAW_GL1PACKET_H

#include "OfflinePacketv1.h"

#include <cstdint>
#include <limits>

class Gl1Packet : public OfflinePacketv1
{
 public:
  Gl1Packet() = default;
  ~Gl1Packet() override = default;
  using OfflinePacket::FillFrom;
  using OfflinePacket::lValue;
  virtual void setScaler(int /*iscal*/, int /*index*/, uint64_t /*lval*/) { return; }
  virtual uint64_t getScaler(int /*iscal*/, int /*index*/) const { return 0; }
  virtual void setGl1pScaler(int /*iscal*/, int /*index*/, uint64_t /*lval*/) { return; }
  virtual uint64_t getGl1pScaler(int /*iscal*/, int /*index*/) const { return 0; }
  virtual void setPacketNumber(const unsigned int /* i */) { return; }
  virtual unsigned int getPacketNumber() const { return 0; }
  virtual void setBunchNumber(const uint64_t /*i*/) { return; }
  virtual uint64_t getBunchNumber() const { return 0; }
  virtual void setTriggerInput(const uint64_t /*i*/) { return; }
  virtual uint64_t getTriggerInput() const { return 0; }
  virtual void setLiveVector(const uint64_t /*i*/) { return; }
  virtual uint64_t getLiveVector() const { return 0; }
  virtual void setTriggerVector(const uint64_t /*i*/) { return; }  // deprecated
  virtual uint64_t getTriggerVector() const { return 0; }          // deprecated
  virtual void setScaledVector(const uint64_t /*i*/) { return; }
  virtual uint64_t getScaledVector() const { return 0; }
  virtual void setGTMBusyVector(const uint64_t /*i*/) { return; }
  virtual uint64_t getGTMBusyVector() const { return 0; }
  virtual void setGTMAllBusyVector(const uint64_t /*i*/) { return; }
  virtual uint64_t getGTMAllBusyVector() const { return 0; }

  virtual void FillFrom(const Gl1Packet * /*pkt*/) { return; }
  // aggregated i/l Values to avoid duplications in the versioned classes
  int iValue(const int i) const override;
  long long lValue(const int i, const std::string &what) const override;
  long long lValue(const int i, const int j) const override;

 private:
  ClassDefOverride(Gl1Packet, 2)
};

#endif
