#include "G4RootHitContainer.h"

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4HitEval.h>

#include <TClonesArray.h>

#include <limits>
#include <ostream>  // for basic_ostream::operator<<, operator<<

static const int NMAX = 100000;

G4RootHitContainer::G4RootHitContainer()
  : SnglHits(new TClonesArray("PHG4HitEval", NMAX))
{
}

G4RootHitContainer::~G4RootHitContainer()
{
  SnglHits->Clear();
  delete SnglHits;
}

void G4RootHitContainer::Reset()
{
  etotal = std::numeric_limits<float>::quiet_NaN();
  leakage = std::numeric_limits<float>::quiet_NaN();
  event = 0;
  SnglHits->Clear();
  if (SnglHits->GetSize() > NMAX)
  {
    SnglHits->Expand(NMAX);
  }
  return;
}

PHG4Hit *
G4RootHitContainer::AddHit(const PHG4Hit *g4hit)
{
  TClonesArray &cl = *SnglHits;
  int nextindex = SnglHits->GetLast() + 1;
  if (nextindex == SnglHits->GetSize())
  {
    SnglHits->Expand(SnglHits->GetSize() + 10000);
  }
  new (cl[nextindex]) PHG4HitEval(g4hit);
  return (static_cast<PHG4Hit *>(cl[nextindex]));
}

void G4RootHitContainer::identify(std::ostream &os) const
{
  os << "Number of Hits: " << SnglHits->GetLast() << std::endl;
  return;
}
