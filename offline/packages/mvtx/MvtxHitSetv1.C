#include "MvtxHitSetv1.h"

  void 
  MvtxHitSetv1::AddHit(const uint16_t col, const uint16_t row)
  {
    hits_.insert(std::make_pair(col, row));
  }

  ConstRange 
  MvtxHitSetv1::GetHits( void )
  {
    return std::make_pair(hits_.begin(), hits_.end());
  }

  ConstRange 
  MvtxHitSetv1::GetHits(const uint16_t col)
  { 
    return std::make_pair(hits_.lower_bound(col), hits_.upper_bound(col));
  }

  
