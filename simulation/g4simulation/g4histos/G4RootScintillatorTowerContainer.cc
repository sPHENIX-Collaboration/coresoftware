#include "G4RootScintillatorTowerContainer.h"

#include "G4RootScintillatorTower.h"

#include <TClonesArray.h>

#include <cmath>                     // for NAN
#include <ostream>                    // for basic_ostream::operator<<, oper...

using namespace std;

static const int NMAX = 1000;

G4RootScintillatorTowerContainer::G4RootScintillatorTowerContainer()
  : idet(-1)
  , etotal(NAN)
  , eion(NAN)
  , leakage(NAN)
  , event(0)
{
  SnglTowers = new TClonesArray("G4RootScintillatorTower", NMAX);
}

G4RootScintillatorTowerContainer::~G4RootScintillatorTowerContainer()
{
  SnglTowers->Clear();
  delete SnglTowers;
}

void G4RootScintillatorTowerContainer::Reset()
{
  etotal = NAN;
  leakage = NAN;
  event = 0;
  SnglTowers->Clear();
  if (SnglTowers->GetSize() > NMAX)
  {
    SnglTowers->Expand(NMAX);
  }
  return;
}

G4RootScintillatorTower *
G4RootScintillatorTowerContainer::AddTower(const RawTower &tower)
{
  TClonesArray &cl = *SnglTowers;
  int nextindex = SnglTowers->GetLast() + 1;
  if (nextindex == SnglTowers->GetSize())
  {
    SnglTowers->Expand(SnglTowers->GetSize() + 10000);
  }
  new (cl[nextindex]) G4RootScintillatorTower(tower);
  return (static_cast<G4RootScintillatorTower *>(cl[nextindex]));
}

void G4RootScintillatorTowerContainer::identify(ostream &os) const
{
  os << "Number of Hits: " << SnglTowers->GetLast() << endl;
  return;
}
