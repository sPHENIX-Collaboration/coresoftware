#ifndef FUN4ALLRAW_MICROMEGASRAWTHIT_H
#define FUN4ALLRAW_MICROMEGASRAWTHIT_H

#include <phool/PHObject.h>

#include <limits>

class MicromegasRawHit : public PHObject
{
 public:
  MicromegasRawHit() = default;
  virtual ~MicromegasRawHit() = default;

  virtual uint64_t get_bco() const { return std::numeric_limits<uint64_t>::max(); }
  virtual void set_bco(const uint64_t) { return; }

  virtual uint64_t get_gtm_bco() const { return std::numeric_limits<uint64_t>::max(); }
  virtual void set_gtm_bco(const uint64_t) { return; }

  virtual int32_t get_packetid() const { return std::numeric_limits<int32_t>::max(); }
  virtual void set_packetid(const int32_t) { return; }

  virtual uint16_t get_fee() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_fee(const uint16_t) { return; }

  virtual uint16_t get_channel() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_channel(const uint16_t) { return; }

  virtual uint16_t get_sampaaddress() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_sampaaddress(const uint16_t) { return; }

  virtual uint16_t get_sampachannel() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_sampachannel(const uint16_t) { return; }

  // index of the first sample with data
  virtual uint16_t get_sample_begin() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_sample_begin(const uint16_t) {}

  // index of the next to last sample with data
  virtual uint16_t get_sample_end() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_sample_end(const uint16_t) {}

  // adc value for a given sample index
  virtual uint16_t get_adc(uint16_t /*sample*/) const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_adc(uint16_t /*sample*/, const uint16_t) { return; }

 private:
  ClassDefOverride(MicromegasRawHit, 1)
};

#endif
