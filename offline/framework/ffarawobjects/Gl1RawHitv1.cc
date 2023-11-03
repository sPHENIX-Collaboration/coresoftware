#include "Gl1RawHitv1.h"

Gl1RawHitv1::Gl1RawHitv1(Gl1RawHit *gl1hit)
{
  set_bco(gl1hit->get_bco());
}

void Gl1RawHitv1::Reset()
{
  bco = std::numeric_limits<uint64_t>::max();
}

void Gl1RawHitv1::identify(std::ostream &os) const
{
  os << "BCO: 0x" << std::hex << bco << std::dec << std::endl;
}
