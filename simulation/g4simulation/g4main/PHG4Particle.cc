#include "PHG4Particle.h"

using namespace std;

void PHG4Particle::identify(ostream& os) const
{
  cout << "calling virtual PHG4Particle base class" << endl;
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
