#include "MvtxHitSetv1.h"


MvtxHitSetv1::MvtxHitSetv1()
  : TrkrHitSetv1()
{

}

void
MvtxHitSetv1::identify(std::ostream& os) const
{
  os << "----- MvtxHitSetv1 -----" << std::endl
     << "       hitid: 0x" << std::hex << GetHitSetKey() << std::dec << std::endl
     << "     truthid: 0x" << std::hex << GetTruthMapKey() << std::dec << std::endl
     << "       nhits: " << hits_.size() << std::endl;

  for ( ConstIterator itr = hits_.begin(); itr != hits_.end(); ++itr)
    os << "            col:" << itr->first << " row:" << itr->second << std::endl;

  os << "------------------------" << std::endl;
}

void
MvtxHitSetv1::Reset()
{
  TrkrHitSetv1::Reset();
  hits_.clear();
}

void
MvtxHitSetv1::print() const
{
  identify(std::cout);
}

void
MvtxHitSetv1::AddHit(const uint16_t col, const uint16_t row)
{
  hits_.insert(std::make_pair(col, row));
}

int
MvtxHitSetv1::RemoveHit(const uint16_t col, const uint16_t row)
{
  for ( Iterator iter = hits_.lower_bound(col);
        iter != hits_.upper_bound(col);
        iter++)
  {
    if ( iter->second == row )
    {
      hits_.erase(iter);
      break;
    }
  }
  return 0;
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


