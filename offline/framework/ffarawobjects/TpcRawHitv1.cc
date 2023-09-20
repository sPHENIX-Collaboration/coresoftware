#include "TpcRawHitv1.h"

TpcRawHitv1::TpcRawHitv1(TpcRawHit *tpchit)
{
  set_bco(tpchit->get_bco());
  set_packetid(tpchit->get_packetid());
  set_fee(tpchit->get_fee());
  set_channel(tpchit->get_channel());
}

void TpcRawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}
