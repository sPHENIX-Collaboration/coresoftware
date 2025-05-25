#include "MicromegasRawHitContainerv2.h"
#include "MicromegasRawHitv2.h"

#include <TClonesArray.h>

static constexpr int NHITS = 100;

MicromegasRawHitContainerv2::MicromegasRawHitContainerv2()
  : MicromegasRawHitsTCArray(new TClonesArray("MicromegasRawHitv2", NHITS))
{
}

MicromegasRawHitContainerv2::~MicromegasRawHitContainerv2()
{
  MicromegasRawHitsTCArray->Clear("C");
  delete MicromegasRawHitsTCArray;
}

void MicromegasRawHitContainerv2::Reset()
{
  MicromegasRawHitsTCArray->Clear("C");
  MicromegasRawHitsTCArray->Expand(NHITS);
}

void MicromegasRawHitContainerv2::identify(std::ostream &os) const
{
  os << "MicromegasRawHitContainerv2" << std::endl;
  os << "containing " << MicromegasRawHitsTCArray->GetEntriesFast() << " Micromegas hits" << std::endl;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
  auto *hit = static_cast<MicromegasRawHit *>(MicromegasRawHitsTCArray->At(0));
  if (hit)
  {
    os << "for beam clock: " << std::hex << hit->get_bco() << std::dec << std::endl;
  }
}

int MicromegasRawHitContainerv2::isValid() const
{
  return MicromegasRawHitsTCArray->GetSize();
}

unsigned int MicromegasRawHitContainerv2::get_nhits()
{
  return MicromegasRawHitsTCArray->GetEntriesFast();
}

MicromegasRawHit *MicromegasRawHitContainerv2::AddHit()
{
  MicromegasRawHit *newhit = new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast() + 1]) MicromegasRawHitv2;
  return newhit;
}

MicromegasRawHit *MicromegasRawHitContainerv2::AddHit(MicromegasRawHit *source)
{
  auto *newhit = new ((*MicromegasRawHitsTCArray)[MicromegasRawHitsTCArray->GetLast() + 1]) MicromegasRawHitv2(source);
  return newhit;
}

MicromegasRawHit *MicromegasRawHitContainerv2::get_hit(unsigned int index)
{
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-static-cast-downcast)
  return static_cast<MicromegasRawHit *>(MicromegasRawHitsTCArray->At(index));
}
