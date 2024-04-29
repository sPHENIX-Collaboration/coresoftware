#include "MicromegasRawHitv1.h"

MicromegasRawHitv1::MicromegasRawHitv1(MicromegasRawHit *source)
{
  set_bco(source->get_bco());
  set_gtm_bco(source->get_gtm_bco());
  set_packetid(source->get_packetid());
  set_fee(source->get_fee());
  set_channel(source->get_channel());
  set_sampaaddress(source->get_sampaaddress());
  set_sampachannel(source->get_sampachannel());
  set_samples(source->get_samples());

  for (size_t i = 0; i < source->get_samples(); ++i)
  {
    set_adc(i, source->get_adc(i));
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
