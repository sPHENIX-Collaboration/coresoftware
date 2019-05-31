#include "G4RootRawTowerContainer.h"

#include "G4RootRawTower.h"

#include <TClonesArray.h>

#include <cmath>
#include <ostream>           // for basic_ostream::operator<<, operator<<, endl

using namespace std;

static const int NMAX = 1000;

G4RootRawTowerContainer::G4RootRawTowerContainer()
  : etotal(NAN)
  , event(0)
{
  SnglG4RootRawTowers = new TClonesArray("G4RootRawTower", NMAX);
}

G4RootRawTowerContainer::~G4RootRawTowerContainer()
{
  SnglG4RootRawTowers->Clear();
  delete SnglG4RootRawTowers;
}

void G4RootRawTowerContainer::Reset()
{
  etotal = NAN;
  event = 0;
  SnglG4RootRawTowers->Clear();
  if (SnglG4RootRawTowers->GetSize() > NMAX)
  {
    SnglG4RootRawTowers->Expand(NMAX);
  }
  return;
}

G4RootRawTower *
G4RootRawTowerContainer::AddG4RootRawTower(const G4RootRawTower &g4tower)
{
  TClonesArray &cl = *SnglG4RootRawTowers;
  int nextindex = SnglG4RootRawTowers->GetLast() + 1;
  if (nextindex == SnglG4RootRawTowers->GetSize())
  {
    SnglG4RootRawTowers->Expand(SnglG4RootRawTowers->GetSize() + 10000);
  }
  new (cl[nextindex]) G4RootRawTower(g4tower);
  return (static_cast<G4RootRawTower *>(cl[nextindex]));
}

void G4RootRawTowerContainer::identify(ostream &os) const
{
  os << "Number of G4RootRawTowers: " << SnglG4RootRawTowers->GetLast() << endl;
  return;
}
