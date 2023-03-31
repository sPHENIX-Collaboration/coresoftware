#include "BbcNorthSouthV1.h"

BbcNorthSouthV1::BbcNorthSouthV1(const Short_t npmt, const float ncharge, const Float_t meantime)
  : nPmt(npmt)
  , nCharge(ncharge)
  , MeanTime(meantime)
{
}


void BbcNorthSouthV1::identify(std::ostream& out) const
{
  out << "identify yourself: I am a BbcNorthSouthV1 object" << std::endl;
}
