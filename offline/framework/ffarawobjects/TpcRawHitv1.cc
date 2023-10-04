#include "TpcRawHitv1.h"

TpcRawHitv1::TpcRawHitv1(TpcRawHit *tpchit)
{
  set_bco(tpchit->get_bco());
  set_gtm_bco(tpchit->get_gtm_bco());
  set_packetid(tpchit->get_packetid());
  set_fee(tpchit->get_fee());
  set_channel(tpchit->get_channel());
  set_sampaaddress(tpchit->get_sampaaddress());
  set_sampachannel(tpchit->get_sampachannel());
  set_samples(tpchit->get_samples());

  for( size_t i = 0; i < tpchit->get_samples(); ++i )
  { set_adc( i, tpchit->get_adc(i) ); }
}

void TpcRawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}
