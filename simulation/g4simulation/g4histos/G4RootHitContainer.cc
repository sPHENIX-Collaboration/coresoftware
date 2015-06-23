#include "G4RootHitContainer.h"

#include <g4main/PHG4Hitv1.h>
#include <TClonesArray.h>

using namespace std;

ClassImp(G4RootHitContainer)

static const int NMAX = 100000;

G4RootHitContainer::G4RootHitContainer():
  etotal(NAN),
  eion(NAN),
  leakage(NAN),
  event(0)
{
 SnglHits = new TClonesArray("PHG4Hitv1",NMAX);
}

G4RootHitContainer::~G4RootHitContainer()
{
  SnglHits->Clear();
  delete SnglHits;
}

void
G4RootHitContainer::Reset()
{
  etotal = NAN;
  leakage = NAN;
  event = 0;
  SnglHits->Clear();
  if (SnglHits->GetSize() > NMAX)
    {
      SnglHits->Expand(NMAX);
    }
  return;
}



PHG4Hit *
G4RootHitContainer::AddHit(const PHG4Hit &g4hit)
{
  TClonesArray &cl = *SnglHits;
  int nextindex = SnglHits->GetLast() + 1;
  if (nextindex == SnglHits->GetSize())
    {
      SnglHits->Expand(SnglHits->GetSize() + 10000);
    }
  new(cl[nextindex]) PHG4Hitv1(g4hit);
  return (static_cast<PHG4Hit *> (cl[nextindex]));
}

void
G4RootHitContainer::identify(ostream& os ) const
{
  os << "Number of Hits: " << SnglHits->GetLast() << endl;
  return;
}
