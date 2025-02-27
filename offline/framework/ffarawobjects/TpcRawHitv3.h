#ifndef FUN4ALLRAW_TPCRAWTHITv3_H
#define FUN4ALLRAW_TPCRAWTHITv3_H

#include "TpcRawHit.h"

#include <phool/PHObject.h>

#include <cassert>
#include <limits>
#include <utility>
#include <vector>

// NOLINTNEXTLINE(hicpp-special-member-functions)
class TpcRawHitv3 : public TpcRawHit
{
 public:
  TpcRawHitv3() = default;
  explicit TpcRawHitv3(TpcRawHit *tpchit);
  TpcRawHitv3(TpcRawHitv3 &&other) noexcept;

  ~TpcRawHitv3() override = default;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  void Clear(Option_t * /*unused*/) override;

  uint64_t get_bco() const override { return bco; }
  // cppcheck-suppress virtualCallInConstructor
  void set_bco(const uint64_t val) override { bco = val; }

  //   uint64_t get_gtm_bco() const override { return gtm_bco; }
  //   // cppcheck-suppress virtualCallInConstructor
  //   void set_gtm_bco(const uint64_t val) override { gtm_bco = val; }

  int32_t get_packetid() const override { return packetid; }
  // cppcheck-suppress virtualCallInConstructor
  void set_packetid(const int32_t val) override { packetid = val; }

  uint16_t get_fee() const override { return fee; }
  // cppcheck-suppress virtualCallInConstructor
  void set_fee(const uint16_t val) override { fee = val; }

  uint16_t get_channel() const override { return channel; }
  // cppcheck-suppress virtualCallInConstructor
  void set_channel(const uint16_t val) override { channel = val; }

  uint16_t get_sampaaddress() const override
  {
    return static_cast<uint16_t>(channel >> 5U) & 0xfU;
  }
  //   // cppcheck-suppress virtualCallInConstructor
  //   void set_sampaaddress(const uint16_t val) override { sampaaddress = val; }

  uint16_t get_sampachannel() const override { return channel & 0x1fU; }
  //   // cppcheck-suppress virtualCallInConstructor
  //   void set_sampachannel(const uint16_t val) override { sampachannel = val; }

  uint16_t get_samples() const override { return 1024U; }
  // cppcheck-suppress virtualCallInConstructor
  //   void set_samples(const uint16_t val) override
  //   {
  //     // assign
  //     samples = val;
  //   }

  uint16_t get_adc(const uint16_t sample) const override;

  // cppcheck-suppress virtualCallInConstructor
  //   void set_adc(const uint16_t sample, const uint16_t val) override
  //   {
  //     adcmap[sample] = val;
  //   }
  void move_adc_waveform(const uint16_t start_time, std::vector<uint16_t> &&adc);

  uint16_t get_type() const override { return type; }
  void set_type(const uint16_t i) override { type = i; }

  // uint16_t get_userword() const override { return userword; }
  // void set_userword(const uint16_t i) override { userword = i; }

  // uint16_t get_checksum() const override { return checksum; }
  // void set_checksum(const uint16_t i) override { checksum = i; }

  // uint16_t get_parity() const override { return data_parity; }
  // void set_parity(const uint16_t i) override { data_parity = i; }

  bool get_checksumerror() const override { return checksumerror; }
  void set_checksumerror(const bool b) override { checksumerror = b; }

  bool get_parityerror() const override { return parityerror; }
  void set_parityerror(const bool b) override { parityerror = b; }

  class AdcIteratorv3 : public AdcIterator
  {
   private:
    const std::vector<std::pair<uint16_t, std::vector<uint16_t> > > &m_adc;
    uint16_t m_waveform_index = 0;
    uint16_t m_adc_position_in_waveform_index = 0;

   public:
    // NOLINTNEXTLINE(hicpp-named-parameter)
    explicit AdcIteratorv3(const std::vector<std::pair<uint16_t, std::vector<uint16_t> > > &adc)
      : m_adc(adc)
    {
    }

    void First() override
    {
      m_waveform_index = 0;
      m_adc_position_in_waveform_index = 0;
    }

    void Next() override
    {
      if (IsDone()) return;
      
      // NOLINTNEXTLINE(bugprone-branch-clone)
      if (m_adc_position_in_waveform_index + 1U < m_adc[m_waveform_index].second.size())
      {
        ++m_adc_position_in_waveform_index;
      }
      else
      { 
        // advance to the next valid waveform
        m_adc_position_in_waveform_index = 0;

        while (true)
        {
          ++m_waveform_index;

          if (m_waveform_index >= m_adc.size())
          {
            break;
          }
          if (not m_adc[m_waveform_index].second.empty())
          {
            break;
          }
        }
      }
    }

    bool IsDone() const override { return m_waveform_index >= m_adc.size(); }

    uint16_t CurrentTimeBin() const override
    {
      if (!IsDone())
      {
        return m_adc[m_waveform_index].first + m_adc_position_in_waveform_index;
      }
      return std::numeric_limits<uint16_t>::max();  // Or throw an exception
    }
    uint16_t CurrentAdc() const override
    {
      if (!IsDone())
      {
        // NOLINTNEXTLINE(bugprone-branch-clone)
        if (m_adc_position_in_waveform_index < m_adc[m_waveform_index].second.size())
        {
          return m_adc[m_waveform_index].second[m_adc_position_in_waveform_index];
        }
        else
        {
          return std::numeric_limits<uint16_t>::max();  // Or throw an exception
        }
      }
      return std::numeric_limits<uint16_t>::max();  // Or throw an exception
    }
  };

  AdcIterator *CreateAdcIterator() const override { return new AdcIteratorv3(m_adcData); }

 private:
  uint64_t bco{std::numeric_limits<uint64_t>::max()};
  int32_t packetid{std::numeric_limits<int32_t>::max()};
  uint16_t fee{std::numeric_limits<uint16_t>::max()};
  uint16_t channel{std::numeric_limits<uint16_t>::max()};
  uint16_t type{std::numeric_limits<uint16_t>::max()};
  // uint16_t checksum{std::numeric_limits<uint16_t>::max()};
  // uint16_t data_parity{std::numeric_limits<uint16_t>::max()};

  bool checksumerror{true};
  bool parityerror{true};

  //! adc waveform std::vector< uint16_t > for each start time uint16_t
  std::vector<std::pair<uint16_t, std::vector<uint16_t> > > m_adcData;

  ClassDefOverride(TpcRawHitv3, 2)
};

#endif
