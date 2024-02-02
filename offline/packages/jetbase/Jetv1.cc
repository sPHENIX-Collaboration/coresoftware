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
std::vector<float> DummyJetPropVecv1;

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

void Jetv1::not_in_v1_msg(const std::string& method_name, std::ostream& os) const {
  os << " warning: Method Jet::"<<method_name <<"() not implemented in Jetv1" << std::endl;
}

std::vector<float>& Jetv1::get_property_vec() {
  not_in_v1_msg("get_property_vec()");
  return DummyJetPropVecv1;
}


// inline float Jetv1::get_prop_by_index(unsigned int /*index*/) const
// {
//   not_in_v1_msg("get_prop_by_index()");
//   return NAN;
// }

//inline void Jetv1::set_prop_by_index(unsigned int /*index*/, float /*value*/)
//{
//  not_in_v1_msg("set_prop_by_index()");
//  return;
//}

void Jetv1::insert_comp(Jet::SRC /**/, unsigned int /**/, bool /**/) 
{
  not_in_v1_msg("insert_comp(src,unsigned int, bool)");
}

void Jetv1::insert_comp(Jet::TYPE_comp_vec&) 
{
  not_in_v1_msg("insert_comp(TYPE_comp_vec&)");
}

void Jetv1::insert_comp(Jet::TYPE_comp_vec&, bool) 
{
  not_in_v1_msg("insert_comp(TYPE_comp_vec&, bool)");
}

void Jetv1::set_comp_sort_flag(bool /**/) 
{
  not_in_v1_msg("set_comp_sort_flag");
}

size_t Jetv1::num_comp(Jet::SRC /**/) 
{
  not_in_v1_msg("num_comp");
  return 0;
}

void Jetv1::print_comp(std::ostream& /**/, bool /**/)
{
  not_in_v1_msg("print_comp");
}

std::vector<Jet::SRC> Jetv1::comp_src_vec() 
{
  not_in_v1_msg("print_comp");
  return {};
}

std::map<Jet::SRC, size_t> Jetv1::comp_src_sizemap()
{
  not_in_v1_msg("comp_src_sizemap");
  return {};
}

Jet::ITER_comp_vec Jetv1::comp_begin(Jet::SRC /**/)
{
  not_in_v1_msg("comp_begin");
  return Jet::comp_begin(Jet::SRC::VOID);
}

Jet::ITER_comp_vec Jetv1::comp_end(Jet::SRC /**/) 
{
  not_in_v1_msg("comp_end");
  return Jet::comp_end(Jet::SRC::VOID);
}

Jet::ITER_comp_vec Jetv1::comp_begin()
{
  not_in_v1_msg("comp_begin");
  return Jet::comp_begin();
}
Jet::ITER_comp_vec Jetv1::comp_end() 
{
  not_in_v1_msg("comp_end");
  return Jet::comp_end();
}
Jet::TYPE_comp_vec& Jetv1::get_comp_vec()
{
  not_in_v1_msg("get_comp_vec");
  return Jet::get_comp_vec();
}
void Jetv1::resize_properties(size_t /**/) 
{
  not_in_v1_msg("resize_properties()");
}
