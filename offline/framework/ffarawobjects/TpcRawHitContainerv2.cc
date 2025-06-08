#include "TpcRawHitContainerv2.h"
#include "TpcRawHitv2.h"

#include <TClonesArray.h>

static const int NTPCHITS = 10000;

TpcRawHitContainerv2::TpcRawHitContainerv2()
  : TpcRawHitsTCArray(new TClonesArray("TpcRawHitv2", NTPCHITS))
{
}

TpcRawHitContainerv2::~TpcRawHitContainerv2()
{
  TpcRawHitsTCArray->Clear("C");
  delete TpcRawHitsTCArray;
}

void TpcRawHitContainerv2::Reset()
{
  TpcRawHitsTCArray->Clear("C");
  TpcRawHitsTCArray->Expand(NTPCHITS);
}

void TpcRawHitContainerv2::identify(std::ostream &os) const
{
  os << "TpcRawHitContainerv2" << std::endl;
  os << "containing " << TpcRawHitsTCArray->GetEntriesFast() << " Tpc hits" << std::endl;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
  TpcRawHit *tpchit = static_cast<TpcRawHit *>(TpcRawHitsTCArray->At(0));
  if (tpchit)
  {
    os << "for beam clock: " << std::hex << tpchit->get_bco() << std::dec << std::endl;
  }
}

int TpcRawHitContainerv2::isValid() const
{
  return TpcRawHitsTCArray->GetSize();
}

unsigned int TpcRawHitContainerv2::get_nhits()
{
  return TpcRawHitsTCArray->GetEntriesFast();
}

TpcRawHit *TpcRawHitContainerv2::AddHit()
{
  TpcRawHit *newhit = new ((*TpcRawHitsTCArray)[TpcRawHitsTCArray->GetLast() + 1]) TpcRawHitv2();
  return newhit;
}

TpcRawHit *TpcRawHitContainerv2::AddHit(TpcRawHit *tpchit)
{
  TpcRawHit *newhit = new ((*TpcRawHitsTCArray)[TpcRawHitsTCArray->GetLast() + 1]) TpcRawHitv2(tpchit);
  return newhit;
}

TpcRawHit *TpcRawHitContainerv2::get_hit(unsigned int index)
{
  return (TpcRawHit *) TpcRawHitsTCArray->At(index);
}
