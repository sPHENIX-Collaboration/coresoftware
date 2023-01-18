#include "BbcPmtHitV1.h"
#include <iostream>

ClassImp(BbcPmtHitV1)

using namespace std;

BbcPmtHitV1::BbcPmtHitV1(const Short_t ipmt, const Float_t iadc, const Float_t itdc0, const Float_t itdc1) :
  pmt(ipmt),
  adc(iadc),
  tdc0(itdc0),
  tdc1(itdc1)
{
}


void BbcPmtHitV1::identify(ostream& out) const
{
  out << "identify yourself: I am a BbcPmtHitV1 object" << endl;
}
