#include "Jetv3.h"

#include <limits>

class PHObject;

Jetv3::Jetv3(unsigned int n_prop)
  : Jetv2(n_prop)
{
}

void Jetv3::identify(std::ostream& os) const
{
  Jetv2::identify(os);
  os << " (Jetv3) corrected pT [GeV]: calib " << _pt_calib
     << ", insitu " << _pt_insitu << " (up " << _pt_insitu_up << ", down " << _pt_insitu_down << ")"
     << ", smear " << _pt_smear << " (up " << _pt_smear_up << ", down " << _pt_smear_down << ")" << std::endl;
}

void Jetv3::Reset()
{
  Jetv2::Reset();
  _pt_calib = std::numeric_limits<float>::quiet_NaN();
  _pt_insitu = std::numeric_limits<float>::quiet_NaN();
  _pt_insitu_up = std::numeric_limits<float>::quiet_NaN();
  _pt_insitu_down = std::numeric_limits<float>::quiet_NaN();
  _pt_smear = std::numeric_limits<float>::quiet_NaN();
  _pt_smear_up = std::numeric_limits<float>::quiet_NaN();
  _pt_smear_down = std::numeric_limits<float>::quiet_NaN();
}

PHObject* Jetv3::CloneMe() const
{
  Jet* jet = new Jetv3(*this);
  return jet;
}
