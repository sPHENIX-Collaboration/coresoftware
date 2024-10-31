#include "MicromegasRawHitv2.h"

MicromegasRawHitv2::MicromegasRawHitv2(MicromegasRawHit *tpchit)
{
  MicromegasRawHitv2::set_bco(tpchit->get_bco());
  MicromegasRawHitv2::set_gtm_bco(tpchit->get_gtm_bco());
  MicromegasRawHitv2::set_packetid(tpchit->get_packetid());
  MicromegasRawHitv2::set_fee(tpchit->get_fee());
  MicromegasRawHitv2::set_channel(tpchit->get_channel());
  MicromegasRawHitv2::set_sampaaddress(tpchit->get_sampaaddress());
  MicromegasRawHitv2::set_sampachannel(tpchit->get_sampachannel());
  MicromegasRawHitv2::set_samples(tpchit->get_samples());

  for (size_t i = 0; i < tpchit->get_samples(); ++i)
  {
    uint16_t adcval = tpchit->get_adc(i);
    if (adcval > 0)
    {
      MicromegasRawHitv2::set_adc(i, tpchit->get_adc(i));
    }
  }
}

void MicromegasRawHitv2::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}

uint16_t MicromegasRawHitv2::get_adc(const uint16_t sample) const
{
  const auto adc = adcmap.find(sample);
  if (adc != adcmap.end())
  {
    return adc->second;
  }
  return 0;
}

void MicromegasRawHitv2::Clear(Option_t * /*unused*/)
{
  adcmap.clear();
}
