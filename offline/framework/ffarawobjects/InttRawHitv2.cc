#include "InttRawHitv2.h"

InttRawHitv2::InttRawHitv2(InttRawHit *intthit)
 : InttRawHitv1(intthit)
{
  InttRawHitv2::set_event_counter(intthit->get_event_counter());
}

void InttRawHitv2::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
  os << "packet id: " << packetid << std::endl;
}
