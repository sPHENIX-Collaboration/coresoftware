#include "TrkrHitSetv1.h"

#include <phool/phool.h>

#include <iostream>

TrkrHitSetv1::TrkrHitSetv1()
  : hitset_key_(TrkrDefs::HITSETKEYMAX)
  , truth_map_key_(UINT64_MAX)
  , hit_clus_map_key_(UINT64_MAX)
{
}

void TrkrHitSetv1::print() const
{
  identify(std::cout);
}

void TrkrHitSetv1::Reset()
{
  hitset_key_ = TrkrDefs::HITSETKEYMAX;
  truth_map_key_ = UINT64_MAX;
  return;
}

void TrkrHitSetv1::identify(std::ostream& os) const
{
  os << "TrkrHitSetv1: " << std::endl
     << "       hitid: 0x" << std::hex << GetHitSetKey() << std::dec << std::endl
     << "     truthid: 0x" << std::hex << GetTruthMapKey() << std::dec << std::endl;
}
