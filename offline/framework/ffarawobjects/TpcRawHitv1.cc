#include "TpcRawHitv1.h"

TpcRawHitv1::TpcRawHitv1(TpcRawHit *tpchit)
{
  TpcRawHitv1::set_bco(tpchit->get_bco());
  TpcRawHitv1::set_gtm_bco(tpchit->get_gtm_bco());
  TpcRawHitv1::set_packetid(tpchit->get_packetid());
  TpcRawHitv1::set_fee(tpchit->get_fee());
  TpcRawHitv1::set_channel(tpchit->get_channel());
  TpcRawHitv1::set_sampaaddress(tpchit->get_sampaaddress());
  TpcRawHitv1::set_sampachannel(tpchit->get_sampachannel());
  TpcRawHitv1::set_samples(tpchit->get_samples());

  for (size_t i = 0; i < tpchit->get_samples(); ++i)
  {
    TpcRawHitv1::set_adc(i, tpchit->get_adc(i));
  }
}

void TpcRawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}

void TpcRawHitv1::Clear(Option_t * /*unused*/)
{
  adc = std::vector<uint16_t>();
}
