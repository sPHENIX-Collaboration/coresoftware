#include "MicromegasRawHitContainerv1.h"
#include "MicromegasRawHitv1.h"

#include <TClonesArray.h>

static constexpr int NHITS = 100;

MicromegasRawHitContainerv1::MicromegasRawHitContainerv1()
{
  MicromegasRawHitsTCArray = new TClonesArray("MicromegasRawHitv1", NHITS);
}

MicromegasRawHitContainerv1::~MicromegasRawHitContainerv1()
{
  delete MicromegasRawHitsTCArray;
}

void MicromegasRawHitContainerv1::Reset()
{
  MicromegasRawHitsTCArray->Clear();
  MicromegasRawHitsTCArray->Expand(NHITS);
}

void MicromegasRawHitContainerv1::identify(std::ostream &os) const
{
  os << "MicromegasRawHitContainerv1" << std::endl;
  os << "containing " << MicromegasRawHitsTCArray->GetEntriesFast() << " Micromegas hits" << std::endl;
  auto hit = static_cast<MicromegasRawHit *>(MicromegasRawHitsTCArray->At(0));
  if (hit)
  {
    os << "for beam clock: " << std::hex << hit->get_bco() << std::dec << std::endl;
  }
}

int MicromegasRawHitContainerv1::isValid() const
{
  return MicromegasRawHitsTCArray->GetSize();
}

unsigned int MicromegasRawHitContainerv1::get_nhits()
{
  return MicromegasRawHitsTCArray->GetEntriesFast();
}

MicromegasRawHit *MicromegasRawHitContainerv1::AddHit()
{
  MicromegasRawHit *newhit = new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast() + 1]) MicromegasRawHitv1();
  return newhit;
}

MicromegasRawHit *MicromegasRawHitContainerv1::AddHit(MicromegasRawHit *source)
{
  auto newhit = new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast() + 1]) MicromegasRawHitv1(source);
  return newhit;
}

MicromegasRawHit *MicromegasRawHitContainerv1::get_hit(unsigned int index)
{
  return static_cast<MicromegasRawHit *>(MicromegasRawHitsTCArray->At(index));
}
