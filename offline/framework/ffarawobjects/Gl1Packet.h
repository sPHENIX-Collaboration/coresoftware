#ifndef FUN4ALLRAW_GL1PACKET_H
#define FUN4ALLRAW_GL1PACKET_H

#include "OfflinePacketv1.h"

#include <limits>

class Gl1Packet : public OfflinePacketv1
{
 public:
  Gl1Packet() = default;
  ~Gl1Packet() override = default;
  using OfflinePacket::FillFrom;
  virtual void setScaler(int /*iscal*/, int /*index*/, uint64_t /*lval*/) { return; }
  virtual void setGl1pScaler(int /*iscal*/, int /*index*/, uint64_t /*lval*/) { return; }
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

  virtual void FillFrom(const Gl1Packet* /*pkt*/) { return; }

 private:
  ClassDefOverride(Gl1Packet, 2)
};

#endif
