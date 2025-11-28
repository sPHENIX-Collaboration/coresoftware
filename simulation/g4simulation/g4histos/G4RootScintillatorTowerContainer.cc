#include "G4RootScintillatorTowerContainer.h"

#include "G4RootScintillatorTower.h"

#include <TClonesArray.h>

#include <limits>
#include <ostream>  // for basic_ostream::operator<<, oper...

static const int NMAX = 1000;

G4RootScintillatorTowerContainer::G4RootScintillatorTowerContainer()
  : SnglTowers(new TClonesArray("G4RootScintillatorTower", NMAX))
{
}

G4RootScintillatorTowerContainer::~G4RootScintillatorTowerContainer()
{
  SnglTowers->Clear();
  delete SnglTowers;
}

void G4RootScintillatorTowerContainer::Reset()
{
  etotal = std::numeric_limits<float>::quiet_NaN();
  leakage = std::numeric_limits<float>::quiet_NaN();
  event = 0;
  SnglTowers->Clear();
  if (SnglTowers->GetSize() > NMAX)
  {
    SnglTowers->Expand(NMAX);
  }
  return;
}

G4RootScintillatorTower *
G4RootScintillatorTowerContainer::AddTower(double towerenergy, int ieta, int iphi)
{
  TClonesArray &cl = *SnglTowers;
  int nextindex = SnglTowers->GetLast() + 1;
  if (nextindex == SnglTowers->GetSize())
  {
    SnglTowers->Expand(SnglTowers->GetSize() + 10000);
  }
  new (cl[nextindex]) G4RootScintillatorTower(towerenergy, ieta, iphi);
  return (static_cast<G4RootScintillatorTower *>(cl[nextindex]));
}

void G4RootScintillatorTowerContainer::identify(std::ostream &os) const
{
  os << "Number of Hits: " << SnglTowers->GetLast() << std::endl;
  return;
}
