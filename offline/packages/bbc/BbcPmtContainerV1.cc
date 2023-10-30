#include "BbcPmtContainerV1.h"
#include "BbcPmtHitV1.h"
#include "BbcReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

static const int NPMTBBCV1 = 128;

BbcPmtContainerV1::BbcPmtContainerV1()
{
  // BbcPmtHit is class for single hit (members: pmt,adc,tdc0,tdc1), do not mix
  // with TClonesArray *BbcPmtHits
  BbcPmtHits = new TClonesArray("BbcPmtHitV1", NPMTBBCV1);
}

BbcPmtContainerV1::~BbcPmtContainerV1()
{
  delete BbcPmtHits;
}

int BbcPmtContainerV1::isValid() const
{
  if (npmt <= 0)
  {
    return 0;
  }
  return 1;
}

void BbcPmtContainerV1::Reset()
{
  BbcPmtHits->Clear();
  npmt = 0;
}

void BbcPmtContainerV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a BbcPmtContainerV1 object" << std::endl;
}
