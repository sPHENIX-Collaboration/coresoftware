#include "TpcRawHitContainerv3.h"
#include "TpcRawHitv3.h"

#include <TClonesArray.h>

static const int NTPCHITS = 10000;

TpcRawHitContainerv3::TpcRawHitContainerv3()
{
  TpcRawHitsTCArray = new TClonesArray("TpcRawHitv3", NTPCHITS);
}

TpcRawHitContainerv3::~TpcRawHitContainerv3()
{
  TpcRawHitsTCArray->Clear("C");
  delete TpcRawHitsTCArray;
}

void TpcRawHitContainerv3::Reset()
{
  TpcRawHitsTCArray->Clear("C");
  TpcRawHitsTCArray->Expand(NTPCHITS);
}

void TpcRawHitContainerv3::identify(std::ostream &os) const
{
  os << "TpcRawHitContainerv3" << std::endl;
  os << "containing " << TpcRawHitsTCArray->GetEntriesFast() << " Tpc hits" << std::endl;
  TpcRawHit *tpchit = static_cast<TpcRawHit *>(TpcRawHitsTCArray->At(0));
  if (tpchit)
  {
    os << "for beam clock: " << std::hex << tpchit->get_bco() << std::dec << std::endl;
  }
}

int TpcRawHitContainerv3::isValid() const
{
  return TpcRawHitsTCArray->GetSize();
}

unsigned int TpcRawHitContainerv3::get_nhits()
{
  return TpcRawHitsTCArray->GetEntriesFast();
}

TpcRawHit *TpcRawHitContainerv3::AddHit()
{
  TpcRawHit *newhit = new ((*TpcRawHitsTCArray)[TpcRawHitsTCArray->GetLast() + 1]) TpcRawHitv3();
  return newhit;
}

TpcRawHit *TpcRawHitContainerv3::AddHit(TpcRawHit *tpchit)
{
  TpcRawHit *newhit = new ((*TpcRawHitsTCArray)[TpcRawHitsTCArray->GetLast() + 1]) TpcRawHitv3(tpchit);
  return newhit;
}

TpcRawHit *TpcRawHitContainerv3::get_hit(unsigned int index)
{
  return (TpcRawHit *) TpcRawHitsTCArray->At(index);
}
