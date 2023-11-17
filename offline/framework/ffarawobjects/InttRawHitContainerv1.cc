#include "InttRawHitContainerv1.h"
#include "InttRawHitv1.h"

#include <TClonesArray.h>

static const int NINTTHITS = 100;

InttRawHitContainerv1::InttRawHitContainerv1()
{
  InttRawHitsTCArray = new TClonesArray("InttRawHitv1",NINTTHITS);
}

InttRawHitContainerv1::~InttRawHitContainerv1()
{
  delete InttRawHitsTCArray;
}


void InttRawHitContainerv1::Reset()
{
 InttRawHitsTCArray->Clear();
 InttRawHitsTCArray->Expand(NINTTHITS);
}

void InttRawHitContainerv1::identify(std::ostream &os) const
{
  os << "InttRawHitContainerv1" << std::endl;
  os << "containing " << InttRawHitsTCArray->GetEntriesFast() << " Intt hits" << std::endl;
  InttRawHit *intthit = static_cast< InttRawHit *> (InttRawHitsTCArray->At(0));
  if (intthit)
   {
     os << "for beam clock: " << std::hex << intthit->get_bco() << std::dec << std::endl;
   }
}

int InttRawHitContainerv1::isValid() const
{
  return InttRawHitsTCArray->GetSize();
}

unsigned int InttRawHitContainerv1::get_nhits()
{
  return InttRawHitsTCArray->GetEntriesFast();
}

InttRawHit *InttRawHitContainerv1::AddHit()
{
  InttRawHit *newhit = new ((*InttRawHitsTCArray)[InttRawHitsTCArray->GetLast()+1]) InttRawHitv1();
  return newhit;
}

InttRawHit *InttRawHitContainerv1::AddHit(InttRawHit *intthit)
{
  InttRawHit *newhit = new ((*InttRawHitsTCArray)[InttRawHitsTCArray->GetLast()+1]) InttRawHitv1(intthit);
  return newhit;
}

InttRawHit *InttRawHitContainerv1::get_hit(unsigned int index)
{
  return (InttRawHit *) InttRawHitsTCArray->At(index);
}
