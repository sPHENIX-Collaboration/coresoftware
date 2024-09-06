#include "TpcRawHitv2.h"

TpcRawHitv2::TpcRawHitv2(TpcRawHit *tpchit)
{
  TpcRawHitv2::set_bco(tpchit->get_bco());
  TpcRawHitv2::set_gtm_bco(tpchit->get_gtm_bco());
  TpcRawHitv2::set_packetid(tpchit->get_packetid());
  TpcRawHitv2::set_fee(tpchit->get_fee());
  TpcRawHitv2::set_channel(tpchit->get_channel());
  TpcRawHitv2::set_sampaaddress(tpchit->get_sampaaddress());
  TpcRawHitv2::set_sampachannel(tpchit->get_sampachannel());
  TpcRawHitv2::set_type(tpchit->get_type());
  TpcRawHitv2::set_userword(tpchit->get_userword());
  TpcRawHitv2::set_checksum(tpchit->get_checksum());
  TpcRawHitv2::set_parity(tpchit->get_parity());
  TpcRawHitv2::set_checksumerror(tpchit->get_checksumerror());
  TpcRawHitv2::set_parityerror(tpchit->get_parityerror());
  TpcRawHitv2::set_samples(tpchit->get_samples());

  for (size_t i = 0; i < tpchit->get_samples(); ++i)
  {
    uint16_t adcval = tpchit->get_adc(i);
    if (adcval > 0)
    {
      TpcRawHitv2::set_adc(i, tpchit->get_adc(i));
    }
  }
}

void TpcRawHitv2::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}

uint16_t TpcRawHitv2::get_adc(const uint16_t sample) const
{
  auto adc = adcmap.find(sample);
  if (adc != adcmap.end())
  {
    return adc->second;
  }
  return 0;
}

void TpcRawHitv2::Clear(Option_t * /*unused*/)
{
  adcmap.clear();
}
