#ifndef FUN4ALLRAW_TPCDIODE_H
#define FUN4ALLRAW_TPCDIODE_H

#include <phool/PHObject.h>

#include <limits>

class TpcDiode : public PHObject
{
 public:
  TpcDiode() = default;
  virtual ~TpcDiode() = default;

  // virtual uint64_t get_bco() const { return std::numeric_limits<uint64_t>::max(); }
  // virtual void set_bco(const uint64_t) { return; }

  // virtual uint64_t get_gtm_bco() const { return std::numeric_limits<uint64_t>::max(); }
  // virtual void set_gtm_bco(const uint64_t) { return; }

  virtual int32_t get_packetid() const { return std::numeric_limits<int32_t>::max(); }
  virtual void set_packetid(const int32_t) { return; }

  // virtual uint16_t get_fee() const { return std::numeric_limits<uint16_t>::max(); }
  // virtual void set_fee(const uint16_t) { return; }

  virtual uint16_t get_channel() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_channel(const uint16_t) { return; }

  // virtual uint16_t get_sampaaddress() const { return std::numeric_limits<uint16_t>::max(); }
  // virtual void set_sampaaddress(const uint16_t) { return; }

  // virtual uint16_t get_sampachannel() const { return std::numeric_limits<uint16_t>::max(); }
  // virtual void set_sampachannel(const uint16_t) { return; }

  virtual uint16_t get_samples() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_samples(const uint16_t) { return; }

  virtual uint16_t get_adc(const uint16_t /*sample*/) const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_adc(const uint16_t /*sample*/, const uint16_t) { return; }

  virtual uint16_t get_maxadc() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_maxadc(const uint16_t) { return; }

  virtual uint16_t get_maxbin() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_maxbin(const uint16_t) { return; }

  virtual double get_integral() const { return std::numeric_limits<double>::max(); }
  virtual void set_integral(const double) { return; }

  virtual uint16_t get_nabovethreshold() const { return std::numeric_limits<uint16_t>::max(); }
  virtual void set_nabovethreshold(const uint16_t) { return; }

  virtual double get_pulsewidth() const { return std::numeric_limits<double>::max(); }
  virtual void set_pulsewidth(const double) { return; }

  // virtual uint16_t get_type() const { return std::numeric_limits<uint16_t>::max(); }
  // virtual void set_type(const uint16_t) { return; }

  // virtual uint16_t get_userword() const { return std::numeric_limits<uint16_t>::max(); }
  // virtual void set_userword(const uint16_t) { return; }

  // virtual uint16_t get_checksum() const { return std::numeric_limits<uint16_t>::max(); }
  // virtual void set_checksum(const uint16_t /*i*/) { return; }

  // virtual uint16_t get_parity() const { return std::numeric_limits<uint16_t>::max(); }
  // virtual void set_parity(const uint16_t /*i*/) { return; }

  // //! FEE waveform CRC check. If true, FEE data transmission from FEE to offline is broken
  // virtual bool get_checksumerror() const { return false; }
  // virtual void set_checksumerror(const bool /*b*/) { return; }

  // //! SAMPA data payload parity check. If true, data from SAMPA is broken, e.g. from SEU
  // virtual bool get_parityerror() const { return false; }
  // virtual void set_parityerror(const bool /*b*/) { return; }

  // class AdcIterator
  // {
  //  public:
  //   AdcIterator() = default;
  //   virtual ~AdcIterator() = default;
  //   virtual void First() = 0;
  //   virtual void Next() = 0;
  //   virtual bool IsDone() const = 0;
  //   virtual uint16_t CurrentTimeBin() const = 0;
  //   virtual uint16_t CurrentAdc() const = 0;
  // };
  // virtual AdcIterator* CreateAdcIterator() const = 0;

 private:
  ClassDefOverride(TpcDiode, 0)
};

#endif
