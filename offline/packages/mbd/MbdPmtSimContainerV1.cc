#include "MbdPmtSimContainerV1.h"
#include "MbdPmtHitV1.h"
#include "MbdReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

static const int NPMTMBDV1 = 128;

MbdPmtSimContainerV1::MbdPmtSimContainerV1()
{
  // MbdPmtHit is class for single hit (members: pmt,adc,tdc0,tdc1), do not mix
  // with TClonesArray *MbdPmtHits
  MbdPmtHits = new TClonesArray("MbdPmtSimHitV1", NPMTMBDV1);
}

MbdPmtSimContainerV1::~MbdPmtSimContainerV1()
{
  delete MbdPmtHits;
}

int MbdPmtSimContainerV1::isValid() const
{
  if (npmt <= 0)
  {
    return 0;
  }
  return 1;
}

void MbdPmtSimContainerV1::Reset()
{
  MbdPmtHits->Clear();
  npmt = 0;
}

void MbdPmtSimContainerV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a MbdPmtSimContainerV1 object" << std::endl;
}
