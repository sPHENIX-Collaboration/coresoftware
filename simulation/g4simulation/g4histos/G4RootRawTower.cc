#include "G4RootRawTower.h"

#include <cmath>
#include <iostream>
#include <limits>

G4RootRawTower::G4RootRawTower(const float ieta, const float iphi, const float e)
  : eta(ieta)
  , phi(iphi)
  , energy(e)
{
}

void G4RootRawTower::Reset()
{
  eta = std::numeric_limits<float>::quiet_NaN();
  phi = std::numeric_limits<float>::quiet_NaN();
  energy = std::numeric_limits<float>::quiet_NaN();
}

int G4RootRawTower::isValid() const
{
  return std::isfinite(get_energy());
}

void G4RootRawTower::identify(std::ostream& os) const
{
  os << "G4RootRawTower: eta: " << eta << ", phi: " << phi
     << " energy=" << get_energy() << std::endl;
}
