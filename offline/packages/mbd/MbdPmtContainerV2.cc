#include "MbdPmtContainerV2.h"
#include "MbdPmtHitV2.h"
#include "MbdReturnCodes.h"

#include <TClonesArray.h>

#include <iostream>

static const int NPMTMBDV2 = 128;

MbdPmtContainerV2::MbdPmtContainerV2() : MbdPmtHits(new TClonesArray("MbdPmtHitV2", NPMTMBDV2))
{
  // MbdPmtHit is class for single hit (members: pmt,adc,tdc0,tdc1), do not mix
  // with TClonesArray *MbdPmtHits
  
}

MbdPmtContainerV2::~MbdPmtContainerV2()
{
  delete MbdPmtHits;
}

int MbdPmtContainerV2::isValid() const
{
  if (npmt <= 0)
  {
    return 0;
  }
  return 1;
}

void MbdPmtContainerV2::Reset()
{
  MbdPmtHits->Clear();
  npmt = 0;
}

void MbdPmtContainerV2::identify(std::ostream &out) const
{
  out << "identify yourself: I am a MbdPmtContainerV2 object" << std::endl;
}
