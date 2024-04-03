#ifndef FUN4ALLRAW_OFFLINEPACKET_H
#define FUN4ALLRAW_OFFLINEPACKET_H

#include <phool/PHObject.h>

#include <iostream>
#include <limits>

class OfflinePacket : public PHObject
{
 public:
  OfflinePacket() = default;
  virtual ~OfflinePacket() = default;
  virtual void FillFrom(const OfflinePacket * /*pkt*/) { return; }

  virtual int getIdentifier() const { return std::numeric_limits<int>::min(); }
  virtual void setIdentifier(const int) { return; }
  virtual int getHitFormat() const { return std::numeric_limits<int>::min(); }
  virtual void setHitFormat(const int) { return; }
  virtual int getEvtSequence() const { return std::numeric_limits<int>::min(); }
  virtual void setEvtSequence(const int) { return; }
  virtual uint64_t getBCO() const { return std::numeric_limits<uint64_t>::max(); }
  virtual void setBCO(const uint64_t) { return; }
  virtual int iValue(const int) const { return std::numeric_limits<int>::min(); }
  virtual int iValue(const int, const std::string &) const { return std::numeric_limits<int>::min(); }
  virtual int iValue(const int, const int) const { return std::numeric_limits<int>::min(); }

  virtual long long lValue(const int /*i*/, const std::string & /*what*/) const { return std::numeric_limits<uint64_t>::max(); }
  virtual long long lValue(const int /*i*/, const int /*j*/) const { return std::numeric_limits<uint64_t>::max(); }

  virtual void dump(std::ostream &os = std::cout) const;

 private:
  ClassDefOverride(OfflinePacket, 1)
};

#endif
