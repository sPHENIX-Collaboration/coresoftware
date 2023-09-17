#include "InttRawHitContainerv1.h"
#include "InttRawHitv1.h"

#include <TClonesArray.h>

static const int NINTTHITS = 100;

InttRawHitContainerv1::InttRawHitContainerv1()
{
  InttRawHitsTCArray = new TClonesArray("InttRawHitsv1",NINTTHITS);
}

InttRawHitContainerv1::~InttRawHitContainerv1()
{
  delete InttRawHitsTCArray;
}


void InttRawHitContainerv1::Reset()
{
InttRawHitsTCArray->Clear();
}

void InttRawHitContainerv1::identify(std::ostream &os) const
{
  os << "InttRawHitContainerv1" << std::endl;
}

int InttRawHitContainerv1::isValid() const
{
  return InttRawHitsTCArray->GetSize();
}

unsigned int InttRawHitContainerv1::get_nhits()
{
  return InttRawHitsTCArray->GetSize();
}

InttRawHit *InttRawHitContainerv1::AddHit()
{
  InttRawHit *newhit = new ((*InttRawHitsTCArray)[InttRawHitsTCArray->GetSize()]) InttRawHitv1();
  return newhit;
}
