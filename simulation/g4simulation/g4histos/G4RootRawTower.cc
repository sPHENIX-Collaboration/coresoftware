#include "G4RootRawTower.h"

#include <cmath>
#include <iostream>

using namespace std;

G4RootRawTower::G4RootRawTower()
  : eta(NAN)
  , phi(NAN)
  , energy(NAN)
{
}

G4RootRawTower::G4RootRawTower(const float ieta, const float iphi, const float e)
  : eta(ieta)
  , phi(iphi)
  , energy(e)
{
}

void G4RootRawTower::Reset()
{
  eta = NAN;
  phi = NAN;
  energy = NAN;
}

int G4RootRawTower::isValid() const
{
  return isfinite(get_energy());
}

void G4RootRawTower::identify(std::ostream& os) const
{
  os << "G4RootRawTower: eta: " << eta << ", phi: " << phi
     << " energy=" << get_energy() << std::endl;
}
