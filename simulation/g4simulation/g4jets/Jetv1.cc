/*!
 * \file Jetv1.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "Jetv1.h"

#include <cmath>
#include <iostream>

using namespace std;

Jetv1::Jetv1()
  : _id(0xFFFFFFFF)
  , _mom()
  , _e(NAN)
  , _comp_ids()
  , _property_map()
{
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
}

void Jetv1::identify(ostream& os) const
{
  os << "---Jet v1-----------------------" << endl;
  os << "jetid: " << get_id() << endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << endl;
  print_property(os);
  for (ConstIter citer = begin_comp(); citer != end_comp(); ++citer)
  {
    cout << citer->first << " -> " << citer->second << endl;
  }
  os << "-----------------------------------------------" << endl;

  return;
}

void Jetv1::Reset()
{
  _id = 0xFFFFFFFF;
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
  _e = NAN;
  _comp_ids.clear();
  _property_map.clear();
}

int Jetv1::isValid() const
{
  if (_id == 0xFFFFFFFF) return 0;
  for (int i = 0; i < 3; ++i)
  {
    if (isnan(_mom[i])) return 0;
  }
  if (isnan(_e)) return 0;
  if (_comp_ids.empty()) return 0;
  return 1;
}

Jet* Jetv1::Clone() const
{
  Jet* jet = new Jetv1(*this);
  return jet;
}

float Jetv1::get_p() const
{
  return sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz());
}

float Jetv1::get_pt() const
{
  return sqrt(get_px() * get_px() + get_py() * get_py());
}

float Jetv1::get_et() const
{
  return get_pt() / get_p() * get_e();
}

float Jetv1::get_eta() const
{
  return asinh(get_pz() / get_pt());
}

float Jetv1::get_phi() const
{
  return atan2(get_py(), get_px());
}

float Jetv1::get_mass() const
{
  // follow CLHEP convention and return negative mass if E^2 - p^2 < 0
  float mass2 = get_mass2();
  if (mass2 < 0)
    return -1 * sqrt(fabs(mass2));
  else
    return sqrt(mass2);
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
    return false;
  else
    return true;
}

float Jetv1::get_property(Jet::PROPERTY prop_id) const
{
  typ_property_map::const_iterator citer = _property_map.find(prop_id);
  if (citer == _property_map.end())
    return NAN;
  else
    return citer->second;
}

void Jetv1::set_property(Jet::PROPERTY prop_id, float value)
{
  _property_map[prop_id] = value;
}

void Jetv1::print_property(ostream& os) const
{
  for (typ_property_map::const_iterator citer = _property_map.begin();
       citer != _property_map.end(); ++citer)
  {
    cout << " ";  //indent

    switch (citer->first)
    {
    case prop_JetCharge:
      cout << "Jet Charge";
      break;
    case prop_BFrac:
      cout << "Jet B-quark fraction";
      break;
    default:
      cout << "Property[" << citer->first << "]";
      break;
    }

    cout << "\t= " << citer->second << endl;
  }
}
