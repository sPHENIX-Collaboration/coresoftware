#include "MicromegasRawHitv1.h"

MicromegasRawHitv1::MicromegasRawHitv1(MicromegasRawHit *source)
{
  MicromegasRawHitv1::set_bco(source->get_bco());
  MicromegasRawHitv1::set_gtm_bco(source->get_gtm_bco());
  MicromegasRawHitv1::set_packetid(source->get_packetid());
  MicromegasRawHitv1::set_fee(source->get_fee());
  MicromegasRawHitv1::set_channel(source->get_channel());
  MicromegasRawHitv1::set_sampaaddress(source->get_sampaaddress());
  MicromegasRawHitv1::set_sampachannel(source->get_sampachannel());
  MicromegasRawHitv1::set_sample_begin(source->get_sample_begin());
  MicromegasRawHitv1::set_sample_end(source->get_sample_end());

  for (size_t i = source->get_sample_begin(); i < source->get_sample_end(); ++i)
  {
    MicromegasRawHitv1::set_adc(i, source->get_adc(i));
  }
}

void MicromegasRawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}

void MicromegasRawHitv1::Clear(Option_t * /*unused*/)
{
  adc = std::vector<uint16_t>();
}
