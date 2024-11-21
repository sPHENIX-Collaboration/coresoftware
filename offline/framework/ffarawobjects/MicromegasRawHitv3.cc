#include "MicromegasRawHitv3.h"

#include <algorithm>
#include <cassert>
#include <iostream>

MicromegasRawHitv3::MicromegasRawHitv3(MicromegasRawHit *source)
{
  static bool once = true;
  if (once)
  {
    once = false;
    std::cout << "MicromegasRawHitv3::MicromegasRawHitv3(MicromegasRawHit *tpchit) - "
              << "WARNING: This moethod is slow and should be avoided as much as possible! Please use the move constructor."
              << std::endl;
  }

  MicromegasRawHitv3::set_bco(source->get_bco());
  MicromegasRawHitv3::set_gtm_bco(source->get_gtm_bco());
  MicromegasRawHitv3::set_packetid(source->get_packetid());
  MicromegasRawHitv3::set_fee(source->get_fee());
  MicromegasRawHitv3::set_channel(source->get_channel());
  MicromegasRawHitv3::set_sampaaddress(source->get_sampaaddress());
  MicromegasRawHitv3::set_sampachannel(source->get_sampachannel());
  MicromegasRawHitv3::set_sample_begin(source->get_sample_begin());
  MicromegasRawHitv3::set_sample_end(source->get_sample_end());

  {
    adc_list_t values;
    for (size_t i = source->get_sample_begin(); i < source->get_sample_end(); ++i)
    { values.push_back( source->get_adc(i)); }

    move_adc_waveform( source->get_sample_begin(), std::move(values) );
  }

}

// cppcheck-suppress accessMoved
MicromegasRawHitv3::MicromegasRawHitv3(MicromegasRawHitv3 &&other) noexcept
  : MicromegasRawHit(other)
  , bco(other.bco)
  , packetid(other.packetid)
  , fee(other.fee)
  , channel(other.channel)
  , m_adcData(std::move(other.m_adcData))
{}

void MicromegasRawHitv3::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}

uint16_t MicromegasRawHitv3::get_adc(const uint16_t sample) const
{
  auto iter = std::find_if( m_adcData.begin(), m_adcData.end(), [sample](const waveform_pair_t& pair)
    { return sample >= pair.first && sample < pair.first+pair.second.size(); } );
  return iter == m_adcData.end() ? 0: iter->second[sample-iter->first];
}

void MicromegasRawHitv3::Clear(Option_t * /*unused*/)
{
  m_adcData.clear();
}

void MicromegasRawHitv3::move_adc_waveform(const uint16_t start_time, std::vector<uint16_t> &&adc)
{
  m_adcData.emplace_back(start_time, adc);
}
