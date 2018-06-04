#include "TrkrHitSet.h"

#include "TrkrHit.h"

#include <phool/phool.h>

#include <iostream>

TrkrHitSet::TrkrHitSet()
  : m_hitSetKey(TrkrDefs::HITSETKEYMAX)
  , m_hits()
{
}

void TrkrHitSet::print() const
{
  identify(std::cout);
}

void TrkrHitSet::Reset()
{
  m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  while ( m_hits.begin() != m_hits.end() )
  {
      delete * m_hits.begin();
      m_hits.erase(m_hits.begin());
  }

  return;
}

void TrkrHitSet::identify(std::ostream& os) const
{
  os << "TrkrHitSetv1: " << std::endl
     << "          id: 0x" << std::hex << getHitSetKey() << std::dec << std::endl
     << "       nhits: " << m_hits.size() << std::endl;

  for ( TrkrHit* hit : m_hits )
  {
      hit->identify(os);
  }
}

unsigned int
TrkrHitSet::addHit(TrkrHit* hit)
{
    m_hits.push_back(hit);
    return m_hits.size() - 1;
}

TrkrHit*
TrkrHitSet::getHit(unsigned int ihit)
{
    // check for input validity
    if ( ihit >= m_hits.size() )
    {
	std::cout << PHWHERE << " - Asked for hit " << ihit
		  << " but only " << size() << " hits available."
		  << " returning nullptr." << std::endl;
	return nullptr;
    }
    return m_hits.at(ihit);
}

TrkrHitSet::ConstRange
TrkrHitSet::getHits()
{
    return std::make_pair(m_hits.begin(), m_hits.end());
}
