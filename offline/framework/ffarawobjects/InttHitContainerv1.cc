#include "InttHitContainerv1.h"
#include "InttHitv1.h"

#include <TClonesArray.h>

static const int NINTTHITS = 100;

InttHitContainerv1::InttHitContainerv1()
{
  InttHitsTCArray = new TClonesArray("InttHitsv1",NINTTHITS);
}

InttHitContainerv1::~InttHitContainerv1()
{
  delete InttHitsTCArray;
}


void InttHitContainerv1::Reset()
{
InttHitsTCArray->Clear();
}

void InttHitContainerv1::identify(std::ostream &os) const
{
  os << "InttHitContainerv1" << std::endl;
}

int InttHitContainerv1::isValid() const
{
  return InttHitsTCArray->GetSize();
}

unsigned int InttHitContainerv1::get_nhits()
{
  return InttHitsTCArray->GetSize();
}

InttHit *InttHitContainerv1::AddHit()
{
  InttHit *newhit = new ((*InttHitsTCArray)[InttHitsTCArray->GetSize()]) InttHitv1();
  return newhit;
}
