#include "G4RootScintillatorTower.h"

// #include <calobase/TowerInfo.h>

#include <iostream>

using namespace std;

G4RootScintillatorTower::G4RootScintillatorTower()
  : row(-1)
  , column(-1)
  , energy(0)
{
}


G4RootScintillatorTower::G4RootScintillatorTower(double towerenergy,int ieta, int iphi)
  : row(iphi)
  , column(ieta)
  , energy(towerenergy)
{

}

void G4RootScintillatorTower::Reset()
{
  row = -1;
  column = -1;
  energy = 0;
}

int G4RootScintillatorTower::isValid() const
{
  return (row >= 0);
}

void G4RootScintillatorTower::identify(std::ostream& os) const
{
  os << "G4RootScintillatorTower: row: " << row << ", column: " << column
     << " energy=" << energy << std::endl;
}
