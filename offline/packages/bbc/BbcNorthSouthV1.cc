#include "BbcNorthSouthV1.h"

BbcNorthSouthV1::BbcNorthSouthV1(const short npmt, const float ncharge, const float meantime)
  : bn(npmt)
  , bq(ncharge)
  , bt(meantime)
{
}

void BbcNorthSouthV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a BbcNorthSouthV1 object" << std::endl;
}
