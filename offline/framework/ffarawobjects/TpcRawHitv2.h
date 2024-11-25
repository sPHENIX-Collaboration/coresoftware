#ifndef FUN4ALLRAW_TPCRAWTHITV2_H
#define FUN4ALLRAW_TPCRAWTHITV2_H

#include "TpcRawHit.h"

#include <phool/PHObject.h>

#include <cassert>
#include <limits>
#include <map>

class TpcRawHitv2 : public TpcRawHit
{
 public:
  TpcRawHitv2() = default;
  TpcRawHitv2(TpcRawHit *tpchit);
  ~TpcRawHitv2() override = default;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  void Clear(Option_t *) override;

  uint64_t get_bco() const override { return bco; }
  // cppcheck-suppress virtualCallInConstructor
  void set_bco(const uint64_t val) override { bco = val; }

  uint64_t get_gtm_bco() const override { return gtm_bco; }
  // cppcheck-suppress virtualCallInConstructor
  void set_gtm_bco(const uint64_t val) override { gtm_bco = val; }

  int32_t get_packetid() const override { return packetid; }
  // cppcheck-suppress virtualCallInConstructor
  void set_packetid(const int32_t val) override { packetid = val; }

  uint16_t get_fee() const override { return fee; }
  // cppcheck-suppress virtualCallInConstructor
  void set_fee(const uint16_t val) override { fee = val; }

  uint16_t get_channel() const override { return channel; }
  // cppcheck-suppress virtualCallInConstructor
  void set_channel(const uint16_t val) override { channel = val; }

  uint16_t get_sampaaddress() const override { return sampaaddress; }
  // cppcheck-suppress virtualCallInConstructor
  void set_sampaaddress(const uint16_t val) override { sampaaddress = val; }

  uint16_t get_sampachannel() const override { return sampachannel; }
  // cppcheck-suppress virtualCallInConstructor
  void set_sampachannel(const uint16_t val) override { sampachannel = val; }

  uint16_t get_samples() const override { return samples; }
  // cppcheck-suppress virtualCallInConstructor
  void set_samples(const uint16_t val) override
  {
    // assign
    samples = val;
  }

  uint16_t get_adc(const uint16_t sample) const override;

  // cppcheck-suppress virtualCallInConstructor
  void set_adc(const uint16_t sample, const uint16_t val) override
  {
    adcmap[sample] = val;
  }

  uint16_t get_type() const override { return type; }
  void set_type(const uint16_t i) override { type = i; }

  uint16_t get_userword() const override { return userword; }
  void set_userword(const uint16_t i) override { userword = i; }

  uint16_t get_checksum() const override { return checksum; }
  void set_checksum(const uint16_t i) override { checksum = i; }

  uint16_t get_parity() const override { return data_parity; }
  void set_parity(const uint16_t i) override { data_parity = i; }

  bool get_checksumerror() const override { return checksumerror; }
  void set_checksumerror(const bool b) override { checksumerror = b; }

  bool get_parityerror() const override { return parityerror; }
  void set_parityerror(const bool b) override { parityerror = b; }

  class AdcIteratorv2 : public AdcIterator
  {
   private:
    const std::map<uint16_t, uint16_t> &m_adc;
    std::map<uint16_t, uint16_t>::const_iterator m_iterator;

   public:
    // NOLINTNEXTLINE(hicpp-member-init)
    explicit AdcIteratorv2(const std::map<uint16_t, uint16_t> &adc)
      : m_adc(adc)
      , m_iterator(adc.begin())
    {
    }

    void First() override { m_iterator = m_adc.begin(); }

    void Next() override { ++m_iterator; }

    bool IsDone() const override { return m_iterator == m_adc.end(); }

    uint16_t CurrentTimeBin() const override
    {
      if (!IsDone())
      {
        return m_iterator->first;
      }
      return std::numeric_limits<uint16_t>::max();  // Or throw an exception
    }
    uint16_t CurrentAdc() const override
    {
      if (!IsDone())
      {
        return m_iterator->second;
      }
      return std::numeric_limits<uint16_t>::max();  // Or throw an exception
    }
  };

  AdcIterator *CreateAdcIterator() const override { return new AdcIteratorv2(adcmap); }

 private:
  uint64_t bco{std::numeric_limits<uint64_t>::max()};
  uint64_t gtm_bco{std::numeric_limits<uint64_t>::max()};
  int32_t packetid{std::numeric_limits<int32_t>::max()};
  uint16_t fee{std::numeric_limits<uint16_t>::max()};
  uint16_t channel{std::numeric_limits<uint16_t>::max()};
  uint16_t sampaaddress{std::numeric_limits<uint16_t>::max()};
  uint16_t sampachannel{std::numeric_limits<uint16_t>::max()};
  uint16_t samples{std::numeric_limits<uint16_t>::max()};
  uint16_t type{std::numeric_limits<uint16_t>::max()};
  uint16_t userword{std::numeric_limits<uint16_t>::max()};
  uint16_t checksum{std::numeric_limits<uint16_t>::max()};
  uint16_t data_parity{std::numeric_limits<uint16_t>::max()};

  bool checksumerror{true};
  bool parityerror{true};

  //! adc value for each sample
  std::map<uint16_t, uint16_t> adcmap;

  ClassDefOverride(TpcRawHitv2, 1)
};

#endif
