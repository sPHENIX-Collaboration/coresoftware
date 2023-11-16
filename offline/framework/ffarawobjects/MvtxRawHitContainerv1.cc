#include "MvtxRawHitContainerv1.h"
#include "MvtxRawHitv1.h"

#include <TClonesArray.h>

static const int NMVTXHITS = 100;

MvtxRawHitContainerv1::MvtxRawHitContainerv1()
{
  MvtxRawHitsTCArray = new TClonesArray("MvtxRawHitv1",NMVTXHITS);
}

MvtxRawHitContainerv1::~MvtxRawHitContainerv1()
{
  delete MvtxRawHitsTCArray;
}


void MvtxRawHitContainerv1::Reset()
{
 MvtxRawHitsTCArray->Clear();
 MvtxRawHitsTCArray->Expand(NMVTXHITS);
}

void MvtxRawHitContainerv1::identify(std::ostream &os) const
{
  os << "MvtxRawHitContainerv1" << std::endl;
  os << "containing " << MvtxRawHitsTCArray->GetEntriesFast() << " Mvtx hits" << std::endl;
  MvtxRawHit *mvtxhit = static_cast< MvtxRawHit *> (MvtxRawHitsTCArray->At(0));
  if (mvtxhit)
   {
     os << "for beam clock: " << std::hex << mvtxhit->get_bco() << std::dec << std::endl;
   }
}

int MvtxRawHitContainerv1::isValid() const
{
  return MvtxRawHitsTCArray->GetSize();
}

unsigned int MvtxRawHitContainerv1::get_nhits()
{
  return MvtxRawHitsTCArray->GetEntriesFast();
}

MvtxRawHit *MvtxRawHitContainerv1::AddHit()
{
  MvtxRawHit *newhit = new ((*MvtxRawHitsTCArray)[MvtxRawHitsTCArray->GetLast()+1]) MvtxRawHitv1();
  return newhit;
}

MvtxRawHit *MvtxRawHitContainerv1::AddHit(MvtxRawHit *mvtxhit)
{
  MvtxRawHit *newhit = new ((*MvtxRawHitsTCArray)[MvtxRawHitsTCArray->GetLast()+1]) MvtxRawHitv1(mvtxhit);
  return newhit;
}

MvtxRawHit *MvtxRawHitContainerv1::get_hit(unsigned int index)
{
  return (MvtxRawHit *) MvtxRawHitsTCArray->At(index);
}
