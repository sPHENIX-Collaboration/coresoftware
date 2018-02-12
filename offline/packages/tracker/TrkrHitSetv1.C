#include "TrkrHitSetv1.h"

#include <phool/phool.h>

#include <iostream>

TrkrHitSetv1::TrkrHitSetv1()
  : hitid(TrkrDefs::HITSETKEYMAX)
  , truthid(UINT64_MAX)
{
}

void TrkrHitSetv1::print() const
{
  identify(std::cout);
}

void TrkrHitSetv1::Reset()
{
  hitid = TrkrDefs::HITSETKEYMAX;
  truthid = UINT64_MAX;
  return;
}

void TrkrHitSetv1::identify(std::ostream& os) const
{
  os << "TrkrHitSetv1: " << std::endl
     << "       hitid: 0x" << std::hex << get_hitid() << std::dec << std::endl
     << "     truthid: 0x" << std::hex << get_truthid() << std::dec << std::endl;
}
