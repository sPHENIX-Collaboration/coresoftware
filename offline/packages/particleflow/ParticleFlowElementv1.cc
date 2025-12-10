#include "ParticleFlowElementv1.h"

#include <cmath>
#include <iostream>
#include <limits>

void ParticleFlowElementv1::identify(std::ostream& os) const
{
  os << "-- ParticleFlowElement v1 : ";
  os << " id: " << get_id() << ", type: " << get_type() << ",";
  os << " (px, py, pz, e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << std::endl;

  return;
}

void ParticleFlowElementv1::Reset()
{
  for (float& i : _mom)
  {
    i = std::numeric_limits<float>::quiet_NaN();
  }
  _e = std::numeric_limits<float>::quiet_NaN();
}

int ParticleFlowElementv1::isValid() const
{
  for (float i : _mom)
  {
    if (std::isnan(i))
    {
      return 0;
    }
  }
  if (std::isnan(_e))
  {
    return 0;
  }

  return 1;
}

float ParticleFlowElementv1::get_p() const
{
  return std::sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz());
}

float ParticleFlowElementv1::get_pt() const
{
  return std::sqrt(get_px() * get_px() + get_py() * get_py());
}

float ParticleFlowElementv1::get_et() const
{
  return get_pt() / get_p() * get_e();
}

float ParticleFlowElementv1::get_eta() const
{
  return std::asinh(get_pz() / get_pt());
}

float ParticleFlowElementv1::get_phi() const
{
  return std::atan2(get_py(), get_px());
}

float ParticleFlowElementv1::get_mass() const
{
  return std::sqrt(get_e() * get_e() - get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz());
}
