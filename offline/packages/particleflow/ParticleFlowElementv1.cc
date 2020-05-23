#include "ParticleFlowElementv1.h"

#include <cmath>
#include <iostream>

ParticleFlowElementv1::ParticleFlowElementv1()
  : _mom() , _e(NAN)
{
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
}

void ParticleFlowElementv1::identify(std::ostream& os) const
{
  os << "-- ParticleFlowElement v1 : ";
  os << " id: " << get_id() << ", ";
  os << " (px, py, pz, e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << std::endl;

  return;
}

void ParticleFlowElementv1::Reset()
{
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
  _e = NAN;
}

int ParticleFlowElementv1::isValid() const
{
  for (int i = 0; i < 3; ++i)
    {
      if (isnan(_mom[i])) return 0;
    }
  if (isnan(_e)) return 0;
  
  return 1;
}

float ParticleFlowElementv1::get_p() const
{
  return sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz());
}

float ParticleFlowElementv1::get_pt() const
{
  return sqrt(get_px() * get_px() + get_py() * get_py());
}

float ParticleFlowElementv1::get_et() const
{
  return get_pt() / get_p() * get_e();
}

float ParticleFlowElementv1::get_eta() const
{
  return asinh(get_pz() / get_pt());
}

float ParticleFlowElementv1::get_phi() const
{
  return atan2(get_py(), get_px());
}

float ParticleFlowElementv1::get_mass() const
{
  
  return sqrt( get_e() * get_e() - get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz() );

}

