#include "MicromegasRawHitv2.h"

MicromegasRawHitv2::MicromegasRawHitv2(MicromegasRawHit *source)
{
  MicromegasRawHitv2::set_bco(source->get_bco());
  MicromegasRawHitv2::set_gtm_bco(source->get_gtm_bco());
  MicromegasRawHitv2::set_packetid(source->get_packetid());
  MicromegasRawHitv2::set_fee(source->get_fee());
  MicromegasRawHitv2::set_channel(source->get_channel());
  MicromegasRawHitv2::set_sampaaddress(source->get_sampaaddress());
  MicromegasRawHitv2::set_sampachannel(source->get_sampachannel());
  MicromegasRawHitv2::set_sample_begin(source->get_sample_begin());
  MicromegasRawHitv2::set_sample_end(source->get_sample_end());

  for (size_t i = source->get_sample_begin(); i < source->get_sample_end(); ++i)
  {
    MicromegasRawHitv2::set_adc(i, source->get_adc(i));
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
