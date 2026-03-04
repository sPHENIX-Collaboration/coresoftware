#include "MbdRawContainerV2.h"
#include "MbdRawHitV2.h"
#include "MbdReturnCodes.h"
#include "MbdDefs.h"

#include <TClonesArray.h>

#include <iostream>

MbdRawContainerV2::MbdRawContainerV2() : MbdRawHits(new TClonesArray("MbdRawHitV2", MbdDefs::MBD_N_PMT))
{
  // MbdRawHit is class for single hit (members: pmt,adc,ttdc,qtdc), do not mix
  // with TClonesArray *MbdRawHits
  
}

MbdRawContainerV2::~MbdRawContainerV2()
{
  delete MbdRawHits;
}

int MbdRawContainerV2::isValid() const
{
  if (npmt <= 0)
  {
    return 0;
  }
  return 1;
}

void MbdRawContainerV2::Reset()
{
  MbdRawHits->Clear();
  npmt = 0;
}

void MbdRawContainerV2::identify(std::ostream &out) const
{
  out << "identify yourself: I am a MbdRawContainerV2 object" << std::endl;
}

//______________________________________
void MbdRawContainerV2::set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk)
{
  evt = ievt;
  clk = iclk;
  femclk = ifemclk;
}

Int_t MbdRawContainerV2::get_evt() const
{
  return evt;
}

UShort_t MbdRawContainerV2::get_clock() const
{
  return clk;
}

UShort_t MbdRawContainerV2::get_femclock() const
{
  return femclk;
}

