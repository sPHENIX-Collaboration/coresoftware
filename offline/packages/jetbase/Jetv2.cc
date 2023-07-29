/*!
 * \file Jetv2.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "Jetv2.h"

#include <phool/phool.h>
#include <algorithm>
#include <cmath>
#include <iostream>

class PHObject;

Jetv2::Jetv2() :
  _sort_ptr { nullptr }
{
  std::fill(std::begin(_mom), std::end(_mom), NAN);
}

Jetv2::Jetv2(unsigned int n_prop)
  : _properties (  n_prop, NAN )
  , _sort_ptr { nullptr }
{
  std::fill(std::begin(_mom), std::end(_mom), NAN);
}

Jetv2::Jetv2(const Jetv2& rhs) 
  : _id               { rhs._id        }
  , _e                { rhs._e         }
  , _is_sorted        { rhs._is_sorted }
  , _sort_ptr         { rhs._sort_ptr  }
{
  std::copy ( rhs._mom, rhs._mom+3, _mom);
  std::copy ( rhs._comp_ids.begin(),   rhs._comp_ids.end(),   _comp_ids.begin() );
  std::copy ( rhs._properties.begin(), rhs._properties.end(), _properties.begin() );
}

void Jetv2::identify(std::ostream& os) const
{
  os << "---Jet v2-----------------------" << std::endl;
  os << "jetid: " << get_id() << std::endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << std::endl;

  os << " Jet Properties:";
  for (auto& val : _properties) { os << " " << val; }
  os << std::endl;

  os << " Jet Components: " << size_comp() << std::endl;;
  return;
}

void Jetv2::print_comp(std::ostream& os, bool single_line) {
  for (auto iter = comp_begin(); iter != comp_end(); ++iter) {
    os << " (" << iter->first << "->" << static_cast<int>(iter->second) << ")";
    if (!single_line) os << std::endl;
  }
  if (single_line) os << std::endl;
}

std::vector<Jet::SRC> Jetv2::comp_src_vec() {
  std::vector<Jet::SRC> vec{};
  if (!_is_sorted) sort_comp_ids();
  auto iter = comp_begin();
  auto iter_end = comp_end();
  if (iter == iter_end) return vec;
  while (iter != iter_end) {
    Jet::SRC src = iter->first;
    vec.push_back(src);
    iter = comp_end(src);
  }
  return vec;
}

std::map<Jet::SRC,size_t> Jetv2::comp_src_sizemap() {
  std::map<Jet::SRC,size_t> sizemap{};
  if (!_is_sorted) sort_comp_ids();
  auto iter = comp_begin();
  auto iter_end = comp_end();
  if (iter == iter_end) return sizemap;
  while (iter != iter_end) {
    Jet::SRC src = iter->first;
    auto iter_ub = comp_end(src);
    sizemap.insert(std::make_pair(src, static_cast<size_t>(iter_ub-iter)));
    iter = iter_ub;
  }
  return sizemap;
}

void Jetv2::Reset()
{
  _id = 0xFFFFFFFF;
  std::fill(std::begin(_mom), std::end(_mom), NAN);
  _e = NAN;
  _comp_ids.clear();
  _properties.clear();
  _is_sorted = false;
}

int Jetv2::isValid() const
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

PHObject* Jetv2::CloneMe() const
{
  Jet* jet = new Jetv2(*this);
  return jet;
}

float Jetv2::get_p() const
{
  return std::sqrt(get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz());
}

float Jetv2::get_pt() const
{
  return std::sqrt(get_px() * get_px() + get_py() * get_py());
}

float Jetv2::get_et() const
{
  return get_pt() / get_p() * get_e();
}

float Jetv2::get_eta() const
{
  return std::asinh(get_pz() / get_pt());
}

float Jetv2::get_phi() const
{
  return std::atan2(get_py(), get_px());
}

float Jetv2::get_mass() const
{
  // follow CLHEP convention and return negative mass if E^2 - p^2 < 0
  float mass2 = get_mass2();
  if (mass2 < 0)
  {
    return -1 * sqrt(std::fabs(mass2));
  }
  return std::sqrt(mass2);
}

float Jetv2::get_mass2() const
{
  float p2 = get_px() * get_px() + get_py() * get_py() + get_pz() * get_pz();
  return get_e() * get_e() - p2;
}

size_t Jetv2::num_comp(Jet::SRC iSRC) {
  if (iSRC==Jet::SRC::VOID) return (comp_end()-comp_begin());
  if (!_is_sorted) { sort_comp_ids(); }
  return (comp_end(iSRC)-comp_begin(iSRC));
}

void Jetv2::insert_comp (SRC iSRC, unsigned int compid)
{ 
  _is_sorted = false; 
  _comp_ids.push_back(std::make_pair(iSRC, compid)); 
}


/* size_t Jetv2::erase_comp(SRC iSRC) { */
/*   size_t n = (comp_end(iSRC) - comp_begin(iSRC)); */
/*   _comp_ids.erase(comp_begin(iSRC), comp_end(iSRC)); */
/*   return n; */
/* } */

void Jetv2::sort_comp_ids() {
    std::sort(_comp_ids.begin(), _comp_ids.end(),
        [](const std::pair<Jet::SRC,unsigned int> &a, 
           const std::pair<Jet::SRC,unsigned int> &b)
        { return a.first < b.first; }
    );
    _is_sorted = true;
}
Jetv2::ITER_comp_vec Jetv2::comp_begin(Jet::SRC iSRC) {
  if (!_is_sorted) sort_comp_ids();
  return std::lower_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

Jetv2::ITER_comp_vec Jetv2::comp_end(Jet::SRC iSRC) {
  if (!_is_sorted) sort_comp_ids();
  return std::upper_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

bool Jetv2::IsEqual(const TObject* obj) const {
  return static_cast<SortFnJetv2*>(_sort_ptr)->IsEqual(
      const_cast<Jetv2*>(this), static_cast<const Jetv2*>(obj));
}

bool SortFnJetv2::IsEqual(const Jetv2* lhs, const Jetv2* rhs) const {
  switch (sort) {
    case Jet::SORT::PT:
      return lhs->get_pt()    == rhs->get_pt();

    case Jet::SORT::E:
      return lhs->get_e()     == rhs->get_e();

    case Jet::SORT::P:
      return lhs->get_p()     == rhs->get_p();

    case Jet::SORT::MASS:
      return lhs->get_mass()  == rhs->get_mass();

    case Jet::SORT::MASS2:
      return lhs->get_mass2() == rhs->get_mass2();

    case Jet::SORT::ETA:
      return lhs->get_eta() == rhs->get_eta();

    case Jet::SORT::PROPERTY:
      return lhs->_properties[index] 
           == rhs->_properties[index];

    default:
      std::cout << PHWHERE << std::endl;
      std::cout << " Unrecognized sorting parameter for Jetv2 types."
        << " Must be Jet::SORT::(PT,E,P,MASS,MASS2) or Jet::SORT::PROPERTY with index within range " << std::endl;
  }
  return true;
}

 int Jetv2::Compare(const TObject* obj) const { 
   return static_cast<SortFnJetv2*>(_sort_ptr)->Compare(
       const_cast<Jetv2*>(this), static_cast<const Jetv2*>(obj)); 
 } 

int SortFnJetv2::Compare(const Jetv2* lhs, const Jetv2* rhs) const {
   switch (sort) {
     case Jet::SORT::PT:
       return intCompare(lhs->get_pt(), rhs->get_pt());
 
     case Jet::SORT::E:
       return intCompare(lhs->get_e(), rhs->get_e());
 
     case Jet::SORT::P:
       return intCompare(lhs->get_p(), rhs->get_p());
 
     case Jet::SORT::MASS:
       return intCompare(lhs->get_mass(), rhs->get_mass());
 
     case Jet::SORT::MASS2:
       return intCompare(lhs->get_mass2(), rhs->get_mass2());
       
     case Jet::SORT::ETA:
       return intCompare(lhs->get_eta(),   rhs->get_eta());
 
     case Jet::SORT::PROPERTY:
       return intCompare(lhs->_properties[index], rhs->_properties[index]);
 
     default:
       std::cout << PHWHERE << std::endl;
       std::cout << " Unrecognized sorting parameter for Jetv2 types."
         << " Must be Jet::SORT::(PT,E,P,MASS,MASS2) or Jet::SORT::PROPERTY with index within range " << std::endl;
   }
   return 0;
 }
