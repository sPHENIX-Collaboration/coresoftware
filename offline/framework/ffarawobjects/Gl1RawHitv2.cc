#include "Gl1RawHitv2.h"

Gl1RawHitv2::Gl1RawHitv2(Gl1RawHit *gl1hit)
{
  set_bco(gl1hit->get_bco());
  setEvtSequence(gl1hit->getEvtSequence());
}

void Gl1RawHitv2::Reset()
{
  Gl1RawHitv1::Reset();
  evtseq = std::numeric_limits<int>::min();
}

void Gl1RawHitv2::identify(std::ostream &os) const
{
  os << "Gl1RawHitv2:" << std::endl;
  os << "Evt Seq: " << getEvtSequence() << std::endl;
  os << "BCO: 0x" << std::hex << get_bco() << std::dec << std::endl;
}
