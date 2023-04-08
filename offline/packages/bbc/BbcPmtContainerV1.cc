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

void BbcPmtContainerV1::AddBbcPmt(const short pmt, const float adc, const float tdc0, const float tdc1)
{
  TClonesArray &Bbchits = *BbcPmtHits;
  new (Bbchits[npmt++]) BbcPmtHitV1(pmt, adc, tdc0, tdc1);
}

short BbcPmtContainerV1::get_pmt(const int iPmt) const
{
  BbcPmtHit *Bbchit = static_cast<BbcPmtHit *> (GetBbcPmtHits()->UncheckedAt(iPmt));
  return ((Bbchit) ? Bbchit->get_pmt() : BbcReturnCodes::BBC_INVALID_SHORT);
}

float BbcPmtContainerV1::get_adc(const int iPmt) const
{
  BbcPmtHit *Bbchit = static_cast<BbcPmtHit *> (GetBbcPmtHits()->UncheckedAt(iPmt));
  return ((Bbchit) ? Bbchit->get_adc() : BbcReturnCodes::BBC_INVALID_FLOAT);
}

float BbcPmtContainerV1::get_tdc0(const int iPmt) const
{
  BbcPmtHit *Bbchit = static_cast<BbcPmtHit *> (GetBbcPmtHits()->UncheckedAt(iPmt));
  return ((Bbchit) ? Bbchit->get_tdc0() : BbcReturnCodes::BBC_INVALID_FLOAT);
}

float BbcPmtContainerV1::get_tdc1(const int iPmt) const
{
  BbcPmtHit *Bbchit = static_cast<BbcPmtHit *> (GetBbcPmtHits()->UncheckedAt(iPmt));
  return ((Bbchit) ? Bbchit->get_tdc1() : BbcReturnCodes::BBC_INVALID_FLOAT);
}

void BbcPmtContainerV1::identify(std::ostream &out) const
{
  out << "identify yourself: I am a BbcPmtContainerV1 object" << std::endl;
}
