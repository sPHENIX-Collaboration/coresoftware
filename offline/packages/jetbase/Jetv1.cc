/*!
 * \file Jetv1.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "Jetv1.h"

#include <algorithm>
#include <cmath>
#include <iostream>

class PHObject;

Jetv1::Jetv1()
{
  std::fill(std::begin(_mom), std::end(_mom), NAN);
}

void Jetv1::identify(std::ostream& os) const
{
  os << "---Jet v1-----------------------" << std::endl;
  os << "jetid: " << get_id() << std::endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << std::endl;
  print_property(os);
  for (ConstIter citer = begin_comp(); citer != end_comp(); ++citer)
  {
    os << citer->first << " -> " << citer->second << std::endl;
  }
  os << "-----------------------------------------------" << std::endl;

  return;
}

void Jetv1::Reset()
{
  _id = 0xFFFFFFFF;
  std::fill(std::begin(_mom), std::end(_mom), NAN);
  _e = NAN;
  _comp_ids.clear();
  _property_map.clear();
}

int Jetv1::isValid() const
{
  if (_id == 0xFFFFFFFF) return 0;
  for (float i : _mom)
  {
    if (std::isnan(i)) return 0;
  }
  if (std::isnan(_e)) return 0;
  if (_comp_ids.empty()) return 0;
  return 1;
}

PHObject* Jetv1::CloneMe() const
{
  Jet* jet = new Jetv1(*this);
  return jet;
}

float Jetv1::get_p() const
{
  return std::sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz());
}

float Jetv1::get_pt() const
{
  return std::sqrt(get_px() * get_px() + get_py() * get_py());
}

float Jetv1::get_et() const
{
  return get_pt() / get_p() * get_e();
}

float Jetv1::get_eta() const
{
  return std::asinh(get_pz() / get_pt());
}

float Jetv1::get_phi() const
{
  return std::atan2(get_py(), get_px());
}

float Jetv1::get_mass() const
{
  // follow CLHEP convention and return negative mass if E^2 - p^2 < 0
  float mass2 = get_mass2();
  if (mass2 < 0)
  {
    return -1 * sqrt(std::fabs(mass2));
  }
  return std::sqrt(mass2);
}

float Jetv1::get_mass2() const
{
  float p2 = get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz();
  return get_e() * get_e() - p2;
}

bool Jetv1::has_property(Jet::PROPERTY prop_id) const
{
  typ_property_map::const_iterator citer = _property_map.find(prop_id);
  if (citer == _property_map.end())
  {
    return false;
  }
  return true;
}

float Jetv1::get_property(Jet::PROPERTY prop_id) const
{
  typ_property_map::const_iterator citer = _property_map.find(prop_id);
  if (citer == _property_map.end())
  {
    return NAN;
  }
  return citer->second;
}

void Jetv1::set_property(Jet::PROPERTY prop_id, float value)
{
  _property_map[prop_id] = value;
}

void Jetv1::print_property(std::ostream& os) const
{
  for (auto citer : _property_map)
  {
    os << " ";  //indent

    switch (citer.first)
    {
    case prop_JetCharge:
      os << "Jet Charge";
      break;
    case prop_BFrac:
      os << "Jet B-quark fraction";
      break;
    default:
      os << "Property[" << citer.first << "]";
      break;
    }

    os << "\t= " << citer.second << std::endl;
  }
}
