#include "TpcRawHitv3.h"

#include <cassert>
#include <iostream>

TpcRawHitv3::TpcRawHitv3(TpcRawHit *tpchit)
{
  static bool once = true;
  if (once)
  {
    once = false;
    std::cout << "TpcRawHitv3::TpcRawHitv3(TpcRawHit *tpchit) - "
              << "WARNING: This moethod is slow and should be avoided as much as possible! Please use the move constructor."
              << std::endl;
  }

  TpcRawHitv3::set_bco(tpchit->get_bco());
  TpcRawHitv3::set_gtm_bco(tpchit->get_gtm_bco());
  TpcRawHitv3::set_packetid(tpchit->get_packetid());
  TpcRawHitv3::set_fee(tpchit->get_fee());
  TpcRawHitv3::set_channel(tpchit->get_channel());
  TpcRawHitv3::set_sampaaddress(tpchit->get_sampaaddress());
  TpcRawHitv3::set_sampachannel(tpchit->get_sampachannel());
  TpcRawHitv3::set_type(tpchit->get_type());
  TpcRawHitv3::set_userword(tpchit->get_userword());
  // TpcRawHitv3::set_checksum(tpchit->get_checksum());
  // TpcRawHitv3::set_parity(tpchit->get_parity());
  TpcRawHitv3::set_checksumerror(tpchit->get_checksumerror());
  TpcRawHitv3::set_parityerror(tpchit->get_parityerror());
  TpcRawHitv3::set_samples(tpchit->get_samples());

  for (size_t i = 0; i < tpchit->get_samples(); ++i)
  {
    uint16_t adcval = tpchit->get_adc(i);
    if (adcval > 0)
    {
      TpcRawHitv3::set_adc(i, tpchit->get_adc(i));
    }
  }
}

// cppcheck-suppress accessMoved
TpcRawHitv3::TpcRawHitv3(TpcRawHitv3 &&other) noexcept
  : TpcRawHit(other)
  , bco(other.bco)
  , packetid(other.packetid)
  , fee(other.fee)
  , channel(other.channel)
  , type(other.type)
  // , userword(other.userword)
  // , checksum(other.checksum)
  // , data_parity(other.data_parity)
  , checksumerror(other.checksumerror)
  , parityerror(other.parityerror)
  , m_adcData(std::move(other.m_adcData))
{
  //   other.bco = std::numeric_limits<uint64_t>::max();
  //   other.packetid = std::numeric_limits<int32_t>::max();
  other.fee = std::numeric_limits<uint16_t>::max();
  other.channel = std::numeric_limits<uint16_t>::max();
  //   other.type = std::numeric_limits<uint16_t>::max();
  //   other.userword = std::numeric_limits<uint16_t>::max();
  //   other.checksum = std::numeric_limits<uint16_t>::max();
  //   other.data_parity = std::numeric_limits<uint16_t>::max();
  other.checksumerror = true;
  other.parityerror = true;
}

void TpcRawHitv3::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << " packet id: " << packetid << std::endl;

  for (const auto &waveform : m_adcData)
  {
    os << " start time: " << waveform.first << " | ADCs: ";
    for (const auto &adc : waveform.second)
    {
      os << adc << " ";
    }
    os << std::endl;
  }
}

uint16_t TpcRawHitv3::get_adc(const uint16_t /*sample*/) const
{
  std::cout << __PRETTY_FUNCTION__
            << " Error: This moethod is slow and should be avoided as much as possible!"
            << std::endl;
  assert(0);
  //   auto adc = adcmap.find(sample);
  //   if (adc != adcmap.end())
  //   {
  //     return adc->second;
  //   }
  return 0;
}

void TpcRawHitv3::Clear(Option_t * /*unused*/)
{
  // quick reset
  fee = std::numeric_limits<uint16_t>::max();
  channel = std::numeric_limits<uint16_t>::max();
  checksumerror = true;
  parityerror = true;

  for (auto &waveform : m_adcData)
  {
    waveform.second.clear();
    waveform.second.shrink_to_fit();
  }
  m_adcData.clear();
  m_adcData.shrink_to_fit();

  // std::cout << __PRETTY_FUNCTION__ << " - m_adcData.capacity = "<<m_adcData.capacity() << std::endl;
}

void TpcRawHitv3::move_adc_waveform(const uint16_t start_time, std::vector<uint16_t> &&adc)
{
  m_adcData.emplace_back(start_time, adc);
}
