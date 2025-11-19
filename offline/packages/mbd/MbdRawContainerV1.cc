#include "MbdRawContainerV1.h"
#include "MbdRawHitV1.h"
#include "MbdReturnCodes.h"
#include "MbdDefs.h"

#include <TClonesArray.h>

#include <iostream>

MbdRawContainerV1::MbdRawContainerV1() : MbdRawHits(new TClonesArray("MbdRawHitV1", MbdDefs::MBD_N_PMT))
{
  // MbdRawHit is class for single hit (members: pmt,adc,ttdc,qtdc), do not mix
  // with TClonesArray *MbdRawHits
  
}

MbdRawContainerV1::~MbdRawContainerV1()
{
  delete MbdRawHits;
}

int MbdRawContainerV1::isValid() const
{
  if (npmt <= 0)
  {
    return 0;
  }
  return 1;
}

void MbdRawContainerV1::Reset()
{
  MbdRawHits->Clear();
  npmt = 0;
}

void MbdRawContainerV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a MbdRawContainerV1 object" << std::endl;
}

//______________________________________
void MbdRawContainerV1::set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk)
{
  evt = ievt;
  clk = iclk;
  femclk = ifemclk;
}

Int_t MbdRawContainerV1::get_evt() const
{
  return evt;
}

UShort_t MbdRawContainerV1::get_clock() const
{
  return clk;
}

UShort_t MbdRawContainerV1::get_femclock() const
{
  return femclk;
}

