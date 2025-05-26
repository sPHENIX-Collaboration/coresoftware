#include "InttRawHitContainerv2.h"
#include "InttRawHitv2.h"

#include <TClonesArray.h>

static const int NINTTHITS = 100;

InttRawHitContainerv2::InttRawHitContainerv2()
  : InttRawHitsTCArray(new TClonesArray("InttRawHitv2", NINTTHITS))
{
}

InttRawHitContainerv2::~InttRawHitContainerv2()
{
  delete InttRawHitsTCArray;
}

void InttRawHitContainerv2::Reset()
{
  InttRawHitsTCArray->Clear();
  InttRawHitsTCArray->Expand(NINTTHITS);
}

void InttRawHitContainerv2::identify(std::ostream &os) const
{
  os << "InttRawHitContainerv2" << std::endl;
  os << "containing " << InttRawHitsTCArray->GetEntriesFast() << " Intt hits" << std::endl;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
  InttRawHit *intthit = static_cast<InttRawHit *>(InttRawHitsTCArray->At(0));
  if (intthit)
  {
    os << "for beam clock: " << std::hex << intthit->get_bco() << std::dec << std::endl;
  }
}

int InttRawHitContainerv2::isValid() const
{
  return InttRawHitsTCArray->GetSize();
}

InttRawHit *InttRawHitContainerv2::AddHit()
{
  InttRawHit *newhit = new ((*InttRawHitsTCArray)[InttRawHitsTCArray->GetLast() + 1]) InttRawHitv2();
  return newhit;
}

InttRawHit *InttRawHitContainerv2::AddHit(InttRawHit *intthit)
{
  InttRawHit *newhit = new ((*InttRawHitsTCArray)[InttRawHitsTCArray->GetLast() + 1]) InttRawHitv2(intthit);
  return newhit;
}

unsigned int InttRawHitContainerv2::get_nhits()
{
  return InttRawHitsTCArray->GetEntriesFast();
}

InttRawHit *InttRawHitContainerv2::get_hit(unsigned int index)
{
  return (InttRawHit *) InttRawHitsTCArray->At(index);
}
