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

Jetv2::Jetv2()
{
  std::fill(std::begin(_mom), std::end(_mom), NAN);
}

Jetv2::Jetv2(unsigned int n_prop) :
  _properties (  n_prop, NAN )
{
  std::fill(std::begin(_mom), std::end(_mom), NAN);
}

Jetv2::Jetv2(const Jetv2& rhs) 
  /* : FnIsEqual  { rhs.FnIsEqual } */
  /* , FnCompare  { rhs.FnCompare } */
  : _id        { rhs._id }
  , _e         { rhs._e }
  , _is_sorted { rhs._is_sorted }
{
  std::copy ( rhs._mom, rhs._mom+3, _mom);
  std::copy ( rhs._comp_ids.begin(),   rhs._comp_ids.end(),   _comp_ids.begin() );
  std::copy ( rhs._properties.begin(), rhs._properties.end(), _properties.begin() );
}

void Jetv2::identify(std::ostream& os) const
{
  os << "---Jet v1-----------------------" << std::endl;
  os << "jetid: " << get_id() << std::endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << std::endl;

  os << " Jet Properties:";
  for (auto& val : _properties) { os << " " << val; }
  os << std::endl;

  for (ConstIter citer = begin_comp(); citer != end_comp(); ++citer)
  {
    os << citer->first << " -> " << citer->second << std::endl;
  }
  os << "-----------------------------------------------" << std::endl;

  return;
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

void Jetv2::print_property(std::ostream& os) const
{
  os << " ::print_property() deprecated in Jetv2; see JetContainer" << std::endl;
  return;
}

size_t Jetv2::cnt_comp(Jet::SRC iSRC) {
  if (!_is_sorted) { sort_comp_ids(); }
  return (comp_end(iSRC)-comp_begin(iSRC));
}


// --deprecated-----------
Jet::typ_comp_ids DummyJetMapv2;

Jet::ConstIter Jetv2::begin_comp() const
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::begin_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::ConstIter Jetv2::lower_bound_comp(Jetv2::SRC /*source*/) const
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::lower_bound_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::ConstIter Jetv2::upper_bound_comp(Jetv2::SRC /*source*/) const
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::upper_bound_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::ConstIter Jetv2::find(Jetv2::SRC /*source*/) const
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::find() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::ConstIter Jetv2::end_comp() const
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::end_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::Iter Jetv2::begin_comp()
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::begin_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::Iter Jetv2::lower_bound_comp(Jetv2::SRC /*source*/)
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::lower_bound_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::Iter Jetv2::upper_bound_comp(Jetv2::SRC /*source*/)
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::upper_bound_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::Iter Jetv2::find(Jetv2::SRC /*source*/)
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::find() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
Jet::Iter Jetv2::end_comp()
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::end_comp() deprecated in Jetv2, getting dummy value" << std::endl;
  return DummyJetMapv2.end();
}
// --deprecated-----------------
size_t Jetv2::count_comp(Jet::SRC /*-*/) const
{
  if (_print_v2_warning) std::cout << PHWHERE 
    << " ::count_comp() deprecated in Jetv2. Use ::cnt_cmp instead" << std::endl;
  return 0;
}

// --deprecated-----------------
void Jetv2::erase_comp(Iter /*first*/, Iter /*last*/) {
  if (_print_v2_warning) {
    std::cout << PHWHERE << " ::erase_comp(Iter,Iter) deprecated in Jetv2" << std::endl;
  }
}
// --deprecated-----------------
void Jetv2::erase_comp(Iter /*first*/) {
  if (_print_v2_warning) { 
    std::cout << PHWHERE << " ::erase_comp(Iter) deprecated in Jetv2" << std::endl;
  }
}
size_t Jetv2::erase_comp(SRC iSRC) {
  if (_print_v2_warning) {
    std::cout << PHWHERE << std::endl;
    std::cout << " Warning: Jetv2 uses vectors for data, and calling"
      << " erase_comp() is very expensive." << std::endl;
  }
  auto ndel = comp_end(iSRC) - comp_begin(iSRC); 
  _comp_ids.erase(comp_begin(iSRC), comp_end(iSRC));
  return ndel;
}

void Jetv2::sort_comp_ids() {
    std::sort(_comp_ids.begin(), _comp_ids.end(),
        [](const std::pair<Jet::SRC,unsigned int> &a, 
           const std::pair<Jet::SRC,unsigned int> &b)
        { return a.first < b.first; }
    );
    _is_sorted = true;
}
Jetv2::TYP_Iter_comp_vec Jetv2::comp_begin(Jet::SRC iSRC) {
  if (!_is_sorted) sort_comp_ids();
  return std::lower_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

Jetv2::TYP_Iter_comp_vec Jetv2::comp_end(Jet::SRC iSRC) {
  if (!_is_sorted) sort_comp_ids();
  return std::upper_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

// -- deprecated 
bool Jetv2::has_property(Jet::PROPERTY /*prop_id*/) const
{
  std::cout << " ::has_property deprecated in Jetv2. Check value in JetContainer " << std::endl;
  return false;
}

float Jetv2::get_property(Jet::PROPERTY index) const
{
  return _properties[index];
}

// -- deprecated 
void Jetv2::set_property(Jet::PROPERTY /**/, float /**/)
{
  std::cout << " ::has_property deprecated in Jetv2. use set_prop_by_index; index per prop from JetContainer " << std::endl;
}


void Jetv2::set_sort_criteria (Jet::SORT which_sort, bool large_to_small, unsigned int index)
{ 
  _which_sort = which_sort; 
  _sort_sign = (large_to_small) ? -1 : 1.;
  _isort_prop_index = index;
};

bool Jetv2::IsEqual(TObject* obj) const {
  Jetv2* rhs { (Jetv2*) obj };
  switch (_which_sort) {
    case Jet::SORT::PT:
      return get_pt()    == rhs->get_pt();

    case Jet::SORT::E:
      return get_e()     == rhs->get_e();

    case Jet::SORT::P:
      return get_p()     == rhs->get_p();

    case Jet::SORT::MASS:
      return get_mass()  == rhs->get_mass();

    case Jet::SORT::MASS2:
      return get_mass2() == rhs->get_mass2();

    case Jet::SORT::ETA:
      return get_eta() == rhs->get_eta();

    case Jet::SORT::PROPERTY:
      return _properties[_isort_prop_index] == rhs->_properties[_isort_prop_index];

    default:
      std::cout << PHWHERE << std::endl;
      std::cout << " Unrecognized sorting parameter for Jetv2 types."
        << " Must be Jet::SORT::(PT,E,P,MASS,MASS2) or Jet::SORT::PROPERTY with index within range " << std::endl;
  }
  return true;
}

int Jetv2::Compare(TObject* obj) {
  Jetv2* rhs { (Jetv2*) obj };
  switch (_which_sort) {
    case Jet::SORT::PT:
      return intCompare(get_pt(), rhs->get_pt());

    case Jet::SORT::E:
      return intCompare(get_e(), rhs->get_e());

    case Jet::SORT::P:
      return intCompare(get_p(), rhs->get_p());

    case Jet::SORT::MASS:
      return intCompare(get_mass(), rhs->get_mass());

    case Jet::SORT::MASS2:
      return intCompare(get_mass2(), rhs->get_mass2());
      
    case Jet::SORT::ETA:
      return intCompare(get_eta(),   rhs->get_eta());

    case Jet::SORT::PROPERTY:
      return intCompare(_properties[_isort_prop_index], rhs->_properties[_isort_prop_index]);

    default:
      std::cout << PHWHERE << std::endl;
      std::cout << " Unrecognized sorting parameter for Jetv2 types."
        << " Must be Jet::SORT::(PT,E,P,MASS,MASS2) or Jet::SORT::PROPERTY with index within range " << std::endl;
  }
  return 0;
}
