#include "PHG4Particlev2.h"

class PHG4Particle;

using namespace std;

PHG4Particlev2::PHG4Particlev2()
  : PHG4Particlev1()
  , trkid(0)
  , vtxid(0)
  , parentid(0)
  , primaryid(0xFFFFFFFF)
  , fe(0.0)
{
}

PHG4Particlev2::PHG4Particlev2(const string &name, const int pid, const double px, const double py, const double pz)
  : PHG4Particlev1(name, pid, px, py, pz)
  , trkid(0)
  , vtxid(0)
  , parentid(0)
  , primaryid(0xFFFFFFFF)
  , fe(0.0)
{
}

PHG4Particlev2::PHG4Particlev2(const PHG4Particle *in)
  : PHG4Particlev1(in)
  , trkid(0)
  , vtxid(0)
  , parentid(0)
  , primaryid(0xFFFFFFFF)
  , fe(0.0)
{
}

void PHG4Particlev2::identify(std::ostream &os) const
{
  if (fname.size() > 0)
  {
    os << "PHG4Particlev2 name: " << fname << ", ";
  }
  else
  {
    os << "PHG4Particlev2 name: missing, ";
  }

  os << "track id: " << trkid
     << ", vtxid: " << vtxid
     << ", parent id: " << parentid
     << ", primary id: " << primaryid
     << ", pid: " << fpid
     << ", px: " << fpx
     << ", py: " << fpy
     << ", pz: " << fpz
     << ", e: " << fe << endl;
  return;
}
