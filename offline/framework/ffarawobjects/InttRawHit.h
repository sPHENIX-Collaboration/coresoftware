#ifndef FUN4ALLRAW_INTTRAWTHIT_H
#define FUN4ALLRAW_INTTRAWTHIT_H

#include <phool/PHObject.h>

#include <limits>

class InttRawHit : public PHObject
{
 public:
  InttRawHit() = default;
  virtual ~InttRawHit() = default;

  virtual uint64_t get_bco() const { return std::numeric_limits<uint64_t>::max(); }
  virtual void set_bco(const uint64_t) { return; }

  virtual int32_t get_packetid() const { return std::numeric_limits<int32_t>::max(); }
  virtual void set_packetid(const int32_t) { return; }

  virtual uint32_t get_word() const { return std::numeric_limits<uint32_t>::max(); }
  virtual void set_word(uint32_t) { return; }

  virtual uint16_t get_fee() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_fee(uint16_t) { return; }

  virtual uint16_t get_channel_id() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_channel_id(uint16_t) { return; }

  virtual uint16_t get_chip_id() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_chip_id(uint16_t) { return; }

  virtual uint16_t get_adc() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_adc(uint16_t) { return; }

  virtual uint16_t get_FPHX_BCO() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_FPHX_BCO(uint16_t) { return; }

  virtual uint16_t get_full_FPHX() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_full_FPHX(uint16_t) { return; }

  virtual uint16_t get_full_ROC() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_full_ROC(uint16_t) { return; }

  virtual uint16_t get_amplitude() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_amplitude(uint16_t) { return; }

 private:
  ClassDefOverride(InttRawHit, 1)
};

#endif
