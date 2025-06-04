#include "PHG4Particlev1.h"

PHG4Particlev1::PHG4Particlev1(const std::string &name, const int pid, const double px, const double py, const double pz)
  : fname(name)
  , fpid(pid)
  , fpx(px)
  , fpy(py)
  , fpz(pz)
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

void PHG4Particlev1::identify(std::ostream &os) const
{
  if (!fname.empty())
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
     << ", barcode: " << barcode << std::endl;
  return;
}
