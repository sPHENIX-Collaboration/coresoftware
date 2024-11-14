#ifndef FUN4ALLRAW_TPCRAWTHITv3_H
#define FUN4ALLRAW_TPCRAWTHITv3_H

#include "MicromegasRawHit.h"

#include <phool/PHObject.h>

#include <cassert>
#include <limits>
#include <utility>
#include <vector>

// NOLINTNEXTLINE(hicpp-special-member-functions)
class MicromegasRawHitv3 : public MicromegasRawHit
{
 public:
  MicromegasRawHitv3() = default;
  explicit MicromegasRawHitv3(MicromegasRawHit*);
  explicit MicromegasRawHitv3(MicromegasRawHitv3 &&other) noexcept;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  void Clear(Option_t */*unused*/) override;

  uint64_t get_bco() const override { return bco; }
  // cppcheck-suppress virtualCallInConstructor
  void set_bco(const uint64_t val) override { bco = val; }

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
  { return static_cast<uint16_t>(channel >> 5U) & 0xfU; }

  uint16_t get_sampachannel() const override { return channel & 0x1fU; }

  // index of the first sample with data
  uint16_t get_sample_begin() const override
  { return m_adcData.empty() ? 0:m_adcData.front().first; }

  // index of the next to last sample with data
  uint16_t get_sample_end() const override
  { return m_adcData.empty() ? 0:m_adcData.back().first+m_adcData.back().second.size(); }

  // get adc value
  uint16_t get_adc(const uint16_t sample) const override;

  // set adc values
  void move_adc_waveform(const uint16_t start_time, std::vector<uint16_t> &&adc);

 private:
  uint64_t bco{std::numeric_limits<uint64_t>::max()};
  int32_t packetid{std::numeric_limits<int32_t>::max()};
  uint16_t fee{std::numeric_limits<uint16_t>::max()};
  uint16_t channel{std::numeric_limits<uint16_t>::max()};
  uint16_t type{std::numeric_limits<uint16_t>::max()};
  uint16_t checksum{std::numeric_limits<uint16_t>::max()};
  uint16_t data_parity{std::numeric_limits<uint16_t>::max()};

  bool checksumerror{true};
  bool parityerror{true};

  //! adc list
  using adc_list_t = std::vector<uint16_t>;

  //! list of adc for a waveform starting at time uint16_t
  using waveform_pair_t = std::pair<uint16_t,adc_list_t>;

  //! list of waveforms
  std::vector<waveform_pair_t> m_adcData;

  ClassDefOverride(MicromegasRawHitv3, 1)
};

#endif
