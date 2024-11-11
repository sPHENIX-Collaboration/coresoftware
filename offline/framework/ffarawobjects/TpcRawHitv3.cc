#include "TpcRawHitv3.h"

#include <iostream>

TpcRawHitv3::TpcRawHitv3(TpcRawHit *tpchit)
{
  static bool once = true;
  if (once)
  {
    once = false;
    std::cout << "TpcRawHitv3::TpcRawHitv3(TpcRawHit *tpchit) - "
              << "WARNING: This moethod is slow and should be avoided as much as possible!"
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
  TpcRawHitv3::set_checksum(tpchit->get_checksum());
  TpcRawHitv3::set_parity(tpchit->get_parity());
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

void TpcRawHitv3::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}

uint16_t TpcRawHitv3::get_adc(const uint16_t sample) const
{
  //   auto adc = adcmap.find(sample);
  //   if (adc != adcmap.end())
  //   {
  //     return adc->second;
  //   }
  return 0;
}

void TpcRawHitv3::Clear(Option_t * /*unused*/)
{
  //   adcmap.clear();
}
