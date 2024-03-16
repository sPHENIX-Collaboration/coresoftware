#include "TpcRawHitContainerv1.h"
#include "TpcRawHitv1.h"

#include <TClonesArray.h>

static const int NTPCHITS = 100;

TpcRawHitContainerv1::TpcRawHitContainerv1()
{
  TpcRawHitsTCArray = new TClonesArray("TpcRawHitv1", NTPCHITS);
}

TpcRawHitContainerv1::~TpcRawHitContainerv1()
{
  TpcRawHitsTCArray->Clear("C");
  delete TpcRawHitsTCArray;
}

void TpcRawHitContainerv1::Reset()
{
  TpcRawHitsTCArray->Clear("C");
  TpcRawHitsTCArray->Expand(NTPCHITS);
}

void TpcRawHitContainerv1::identify(std::ostream &os) const
{
  os << "TpcRawHitContainerv1" << std::endl;
  os << "containing " << TpcRawHitsTCArray->GetEntriesFast() << " Tpc hits" << std::endl;
  TpcRawHit *tpchit = static_cast<TpcRawHit *>(TpcRawHitsTCArray->At(0));
  if (tpchit)
  {
    os << "for beam clock: " << std::hex << tpchit->get_bco() << std::dec << std::endl;
  }
}

int TpcRawHitContainerv1::isValid() const
{
  return TpcRawHitsTCArray->GetSize();
}

unsigned int TpcRawHitContainerv1::get_nhits()
{
  return TpcRawHitsTCArray->GetEntriesFast();
}

TpcRawHit *TpcRawHitContainerv1::AddHit()
{
  TpcRawHit *newhit = new ((*TpcRawHitsTCArray)[TpcRawHitsTCArray->GetLast() + 1]) TpcRawHitv1();
  return newhit;
}

TpcRawHit *TpcRawHitContainerv1::AddHit(TpcRawHit *tpchit)
{
  TpcRawHit *newhit = new ((*TpcRawHitsTCArray)[TpcRawHitsTCArray->GetLast() + 1]) TpcRawHitv1(tpchit);
  return newhit;
}

TpcRawHit *TpcRawHitContainerv1::get_hit(unsigned int index)
{
  return (TpcRawHit *) TpcRawHitsTCArray->At(index);
}
