#include "BbcPmtHit.h"
#include <iostream>

ClassImp(BbcPmtHit)

using namespace std;

void BbcPmtHit::identify(ostream& out) const
{
  out << "identify yourself: I am a BbcPmtHit object" << endl;
}
