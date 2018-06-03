#include "TrkrHitSetv1.h"

#include <phool/phool.h>

#include <iostream>

TrkrHitSetv1::TrkrHitSetv1()
  : m_hitSetKey(TrkrDefs::HITSETKEYMAX)
  , m_truthMapKey(UINT64_MAX)
  , m_hitClusMapKey(UINT64_MAX)
{
}

void TrkrHitSetv1::print() const
{
  identify(std::cout);
}

void TrkrHitSetv1::Reset()
{
  m_hitSetKey = TrkrDefs::HITSETKEYMAX;
  m_truthMapKey = UINT64_MAX;
  return;
}

void TrkrHitSetv1::identify(std::ostream& os) const
{
  os << "TrkrHitSetv1: " << std::endl
     << "       hitid: 0x" << std::hex << getHitSetKey() << std::dec << std::endl
     << "     truthid: 0x" << std::hex << getTruthMapKey() << std::dec << std::endl;
}
