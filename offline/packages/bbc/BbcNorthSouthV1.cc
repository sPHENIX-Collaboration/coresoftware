#include "BbcNorthSouthV1.h"

using namespace std;

ClassImp(BbcNorthSouthV1)

BbcNorthSouthV1::BbcNorthSouthV1(const Short_t npmt, const Float_t ncharge, const Float_t meantime)
{
  nPmt     = npmt;
  nCharge  = ncharge;
  MeanTime = meantime;
}


void BbcNorthSouthV1::identify(ostream& out) const
{
  out << "identify yourself: I am a BbcNorthSouthV1 object" << endl;
}
