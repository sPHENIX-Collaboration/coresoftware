#include "TpcRawHitv2.h"

TpcRawHitv2::TpcRawHitv2(TpcRawHit *tpchit)
{
  set_bco(tpchit->get_bco());
  set_gtm_bco(tpchit->get_gtm_bco());
  set_packetid(tpchit->get_packetid());
  set_fee(tpchit->get_fee());
  set_channel(tpchit->get_channel());
  set_sampaaddress(tpchit->get_sampaaddress());
  set_sampachannel(tpchit->get_sampachannel());
  set_samples(tpchit->get_samples());

  for (size_t i = 0; i < tpchit->get_samples(); ++i)
  {
    uint16_t adcval = tpchit->get_adc(i);
    if (adcval > 0)
    {
      set_adc(i, tpchit->get_adc(i));
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
