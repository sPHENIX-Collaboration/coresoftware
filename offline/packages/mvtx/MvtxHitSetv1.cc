#include "MvtxHitSetv1.h"


MvtxHitSetv1::MvtxHitSetv1()
{

}

void 
MvtxHitSetv1::identify(std::ostream& os) const
{

}

void 
MvtxHitSetv1::Reset()
{

}

void
MvtxHitSetv1::print() const
{

}

void 
MvtxHitSetv1::addHit(const uint16_t col, const uint16_t row)
{
  m_hits.insert(std::make_pair(col, row));
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

  
