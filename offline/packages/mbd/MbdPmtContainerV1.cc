#include "MbdPmtContainerV1.h"
#include "MbdPmtHitV1.h"
#include "MbdReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

static const int NPMTMBDV1 = 128;

MbdPmtContainerV1::MbdPmtContainerV1()
{
  // MbdPmtHit is class for single hit (members: pmt,adc,tdc0,tdc1), do not mix
  // with TClonesArray *MbdPmtHits
  MbdPmtHits = new TClonesArray("MbdPmtHitV1", NPMTMBDV1);
}

MbdPmtContainerV1::~MbdPmtContainerV1()
{
  delete MbdPmtHits;
}

int MbdPmtContainerV1::isValid() const
{
  if (npmt <= 0)
  {
    return 0;
  }
  return 1;
}

void MbdPmtContainerV1::Reset()
{
  MbdPmtHits->Clear();
  npmt = 0;
}

void MbdPmtContainerV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a MbdPmtContainerV1 object" << std::endl;
}
