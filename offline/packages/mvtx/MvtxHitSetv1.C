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
MvtxHitSetv1::AddHit(const uint16_t col, const uint16_t row)
{
  hits_.insert(std::make_pair(col, row));
}

MvtxHitSetv1::ConstRange 
MvtxHitSetv1::GetHits( void )
{
  return std::make_pair(hits_.begin(), hits_.end());
}

MvtxHitSetv1::ConstRange 
MvtxHitSetv1::GetHits(const uint16_t col)
{ 
  return std::make_pair(hits_.lower_bound(col), hits_.upper_bound(col));
}

  
