#include "BbcPmtContainerV1.h"
#include "BbcPmtHitV1.h"
#include "BbcReturnCodes.h"
#include <TClonesArray.h>
#include <iostream>

using namespace std;

static const int NPMTBBCV1 = 128;

ClassImp(BbcPmtContainerV1)

BbcPmtContainerV1::BbcPmtContainerV1()
{
  // BbcPmtHit is class for single hit (members: pmt,adc,tdc0,tdc1), do not mix
  // with TClonesArray *BbcPmtHits
  BbcPmtHits = new TClonesArray("BbcPmtHit", NPMTBBCV1);
  npmt = 0;
}

BbcPmtContainerV1::~BbcPmtContainerV1()
{
  if (BbcPmtHits)
    {
      delete BbcPmtHits;
    }
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

void BbcPmtContainerV1::AddBbcPmt(const Short_t pmt, const Float_t adc, const Float_t tdc0, const Float_t tdc1)
{
  TClonesArray &Bbchits = *BbcPmtHits;
  new(Bbchits[pmt]) BbcPmtHitV1(pmt, adc, tdc0, tdc1);
}

Float_t BbcPmtContainerV1::get_adc(const int iPmt) const
{
  BbcPmtHit *Bbchit = (BbcPmtHit*) GetBbcPmtHits()->UncheckedAt(iPmt);
  return ((Bbchit) ? Bbchit->get_adc() : BBC_INVALID_SHORT);
}

Float_t BbcPmtContainerV1::get_tdc0(const int iPmt) const
{
  BbcPmtHit *Bbchit = (BbcPmtHit*) GetBbcPmtHits()->UncheckedAt(iPmt);
  return ((Bbchit) ? Bbchit->get_tdc0() : BBC_INVALID_SHORT);
}

Float_t BbcPmtContainerV1::get_tdc1(const int iPmt) const
{
  BbcPmtHit *Bbchit = (BbcPmtHit*) GetBbcPmtHits()->UncheckedAt(iPmt);
  return ((Bbchit) ? Bbchit->get_tdc1() : BBC_INVALID_SHORT);
}

void BbcPmtContainerV1::identify(ostream& out) const
{
  out << "identify yourself: I am a BbcPmtContainerV1 object" << endl;
}
