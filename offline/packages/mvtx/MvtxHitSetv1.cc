#include "MvtxHitSetv1.h"


MvtxHitSetv1::MvtxHitSetv1()
  : TrkrHitSetv1()
{

}

void
MvtxHitSetv1::identify(std::ostream& os) const
{
  os << "----- MvtxHitSetv1 -----" << std::endl
     << "       hitid: 0x" << std::hex << getHitSetKey() << std::dec << std::endl
     << "     truthid: 0x" << std::hex << getTruthMapKey() << std::dec << std::endl
     << "       nhits: " << m_hits.size() << std::endl;

  for ( ConstIterator itr = m_hits.begin(); itr != m_hits.end(); ++itr)
    os << "            col:" << itr->first << " row:" << itr->second << std::endl;

  os << "------------------------" << std::endl;
}

void
MvtxHitSetv1::Reset()
{
  TrkrHitSetv1::Reset();
  m_hits.clear();
}

void
MvtxHitSetv1::print() const
{
  identify(std::cout);
}

void
MvtxHitSetv1::addHit(const uint16_t col, const uint16_t row)
{
  m_hits.insert(std::make_pair(col, row));
}

int
MvtxHitSetv1::removeHit(const uint16_t col, const uint16_t row)
{
  for ( Iterator iter = m_hits.lower_bound(col);
        iter != m_hits.upper_bound(col);
        iter++)
  {
    if ( iter->second == row )
    {
      m_hits.erase(iter);
      break;
    }
  }
  return 0;
}

MvtxHitSetv1::ConstRange
MvtxHitSetv1::getHits( void )
{
  return std::make_pair(m_hits.begin(), m_hits.end());
}

MvtxHitSetv1::ConstRange
MvtxHitSetv1::getHits(const uint16_t col)
{
  return std::make_pair(m_hits.lower_bound(col), m_hits.upper_bound(col));
}


