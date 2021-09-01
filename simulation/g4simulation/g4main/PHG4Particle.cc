#include "PHG4Particle.h"

void PHG4Particle::identify(std::ostream& os) const
{
  os << "calling virtual PHG4Particle base class" << std::endl;
  return;
}

bool PHG4Particle::operator==(const PHG4Particle& p) const
{
  if (get_pid() == p.get_pid() &&
      get_px() == p.get_px() &&
      get_py() == p.get_py() &&
      get_pz() == p.get_pz())
  {
    return true;
  }
  return false;
}
