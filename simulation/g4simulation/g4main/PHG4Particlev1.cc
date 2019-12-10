#include "PHG4Particlev1.h"

using namespace std;

PHG4Particlev1::PHG4Particlev1()
  : fname("")
  , fpid(0)
  , fpx(0)
  , fpy(0)
  , fpz(0)
  , barcode(-1)
{
}

PHG4Particlev1::PHG4Particlev1(const string &name, const int pid, const double px, const double py, const double pz)
  : fname(name)
  , fpid(pid)
  , fpx(px)
  , fpy(py)
  , fpz(pz)
  , barcode(-1)
{
}

PHG4Particlev1::PHG4Particlev1(const PHG4Particle *in)
  : fname(in->get_name())
  , fpid(in->get_pid())
  , fpx(in->get_px())
  , fpy(in->get_py())
  , fpz(in->get_pz())
  , barcode(in->get_barcode())
{
}

void PHG4Particlev1::identify(ostream &os) const
{
  if (fname.size() > 0)
  {
    os << "PHG4Particlev1 name: " << fname;
  }
  else
  {
    os << "PHG4Particlev1 name: missing ";
  }
  os << ", pid: " << fpid
     << ", px: " << fpx
     << ", py: " << fpy
     << ", pz: " << fpz
     << ", barcode: " << barcode << endl;
  return;
}
