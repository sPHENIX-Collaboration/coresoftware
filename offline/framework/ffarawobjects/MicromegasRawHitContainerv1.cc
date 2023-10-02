#include "MicromegasRawHitContainerv1.h"
#include "MicromegasRawHitv1.h"

#include <TClonesArray.h>

static const int NTPCHITS = 100;

MicromegasRawHitContainerv1::MicromegasRawHitContainerv1()
{
  MicromegasRawHitsTCArray = new TClonesArray("MicromegasRawHitv1",NTPCHITS);
}

MicromegasRawHitContainerv1::~MicromegasRawHitContainerv1()
{
  delete MicromegasRawHitsTCArray;
}


void MicromegasRawHitContainerv1::Reset()
{
 MicromegasRawHitsTCArray->Clear();
 MicromegasRawHitsTCArray->Expand(NTPCHITS);
}

void MicromegasRawHitContainerv1::identify(std::ostream &os) const
{
  os << "MicromegasRawHitContainerv1" << std::endl;
  os << "containing " << MicromegasRawHitsTCArray->GetEntries() << " Tpc hits" << std::endl;
  MicromegasRawHit *tpchit = static_cast< MicromegasRawHit *> (MicromegasRawHitsTCArray->At(0));
  if (tpchit)
   {
     os << "for beam clock: " << std::hex << tpchit->get_bco() << std::dec << std::endl;
   }
}

int MicromegasRawHitContainerv1::isValid() const
{
  return MicromegasRawHitsTCArray->GetSize();
}

unsigned int MicromegasRawHitContainerv1::get_nhits()
{
  return MicromegasRawHitsTCArray->GetEntries();
}

MicromegasRawHit *MicromegasRawHitContainerv1::AddHit()
{
  MicromegasRawHit *newhit = new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast()+1]) MicromegasRawHitv1();
  return newhit;
}

MicromegasRawHit *MicromegasRawHitContainerv1::AddHit(MicromegasRawHit *tpchit)
{
  MicromegasRawHit *newhit = new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast()+1]) MicromegasRawHitv1(tpchit);
  return newhit;
}

MicromegasRawHit *MicromegasRawHitContainerv1::get_hit(unsigned int index)
{
  return (MicromegasRawHit *) MicromegasRawHitsTCArray->At(index);
}
