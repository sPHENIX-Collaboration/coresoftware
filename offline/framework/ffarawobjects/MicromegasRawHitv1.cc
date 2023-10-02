#include "MicromegasRawHitv1.h"

MicromegasRawHitv1::MicromegasRawHitv1(MicromegasRawHit *tpchit)
{
  set_bco(tpchit->get_bco());
  set_packetid(tpchit->get_packetid());
  set_fee(tpchit->get_fee());
  set_channel(tpchit->get_channel());
  set_sampaaddress(tpchit->get_sampaaddress());
  set_sampachannel(tpchit->get_sampachannel());
  set_samples(tpchit->get_samples());
}

void MicromegasRawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}
