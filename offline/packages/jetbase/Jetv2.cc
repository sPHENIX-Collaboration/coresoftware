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

Jetv2::Jetv2(unsigned int n_prop)
  : _properties(n_prop, NAN)
{
  std::fill(std::begin(_mom), std::end(_mom), NAN);
}

Jetv2::Jetv2(const Jetv2& rhs)
  : _id{rhs._id}
  , _e{rhs._e}
{
  std::copy(rhs._mom, rhs._mom + 3, _mom);
  std::copy(rhs._comp_ids.begin(), rhs._comp_ids.end(), _comp_ids.begin());
  std::copy(rhs._properties.begin(), rhs._properties.end(), _properties.begin());
}

void Jetv2::identify(std::ostream& os) const
{
  os << "---Jet v2-----------------------" << std::endl;
  os << "jetid: " << get_id() << std::endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << std::endl;

  os << " Jet Properties:";
  for (auto& val : _properties)
  {
    os << " " << val;
  }
  os << std::endl;

  os << " Jet Components: " << _comp_ids.size() << std::endl;
  ;
  return;
}

void Jetv2::print_comp(std::ostream& os, bool single_line)
{
  for (auto iter = comp_begin(); iter != comp_end(); ++iter)
  {
    os << " (" << iter->first << "->" << static_cast<int>(iter->second) << ")";
    if (!single_line) os << std::endl;
  }
  if (single_line) os << std::endl;
}

std::vector<Jet::SRC> Jetv2::comp_src_vec()
{
  std::vector<Jet::SRC> vec{};
  auto iter = comp_begin();
  auto iter_end = comp_end();
  if (iter == iter_end) return vec;
  while (iter != iter_end)
  {
    Jet::SRC src = iter->first;
    vec.push_back(src);
    iter = comp_end(src);
  }
  return vec;
}

std::map<Jet::SRC, size_t> Jetv2::comp_src_sizemap()
{
  std::map<Jet::SRC, size_t> sizemap{};
  auto iter = comp_begin();
  auto iter_end = comp_end();
  if (iter == iter_end) return sizemap;
  while (iter != iter_end)
  {
    Jet::SRC src = iter->first;
    auto iter_ub = comp_end(src);
    sizemap.insert(std::make_pair(src, static_cast<size_t>(iter_ub - iter)));
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

size_t Jetv2::num_comp(Jet::SRC iSRC)
{
  if (iSRC == Jet::SRC::VOID) return (comp_end() - comp_begin());
  return (comp_end(iSRC) - comp_begin(iSRC));
}

void Jetv2::insert_comp(SRC iSRC, unsigned int compid)
{
  _comp_ids.push_back(std::make_pair(iSRC, compid));
}

void Jetv2::CompareSRC::sort_comp_ids(Jetv2* jet)
{
  std::sort(jet->_comp_ids.begin(), jet->_comp_ids.end(),
            [](const std::pair<Jet::SRC, unsigned int>& a,
               const std::pair<Jet::SRC, unsigned int>& b)
            { 
            if (a.first == b.first) return a.second < b.second;
            else return a.first < b.first; });
}

Jetv2::ITER_comp_vec Jetv2::comp_begin(Jet::SRC iSRC)
{
  return std::lower_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

Jetv2::ITER_comp_vec Jetv2::comp_end(Jet::SRC iSRC) 
{
  return std::upper_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

void Jetv2::msg_dep_fn(const std::string&  fn_name) const {
  std::cout << " warning: Method Jet::"<<fn_name <<"() deprecated in Jetv2" << std::endl;
}

bool Jetv2::has_property(Jet::PROPERTY /*prop_id*/) const {
    msg_dep_fn("has_property"); 
    return false;
}

float Jetv2::get_property(Jet::PROPERTY /*prop_id*/) const {
    msg_dep_fn("get_property"); 
    return NAN;
}

void Jetv2::set_property(Jet::PROPERTY /**/, float /**/) {
    msg_dep_fn("set_property"); 
}

void Jetv2::print_property(std::ostream& /**/) const {
    msg_dep_fn("print_property"); 
}

bool Jetv2::empty_comp() const {
    msg_dep_fn("empty_comp"); 
    return true;
}

size_t Jetv2::count_comp(Jet::SRC /**/) const {
  msg_dep_fn("count_comp"); 
  return 0;
}

Jet::ConstIter Jetv2::begin_comp() const
{
  msg_dep_fn("begin_comp");
  return Jet::begin_comp();
}
Jet::ConstIter Jetv2::lower_bound_comp(Jet::SRC /**/) const
{
  msg_dep_fn("lower_bound_comp");
  return Jet::lower_bound_comp(Jet::SRC::VOID);
}
Jet::ConstIter Jetv2::upper_bound_comp(Jet::SRC /**/) const
{
  msg_dep_fn("upper_bound_comp");
  return Jet::upper_bound_comp(Jet::SRC::VOID);
}
Jet::ConstIter Jetv2::find(Jet::SRC /**/) const
{
  msg_dep_fn("find");
  return Jet::find(Jet::SRC::VOID);
}
Jet::ConstIter Jetv2::end_comp() const
{
  msg_dep_fn("end_comp");
  return Jet::end_comp();
}
Jet::Iter Jetv2::begin_comp()
{
  msg_dep_fn("begin_comp");
  return Jet::begin_comp();
}
Jet::Iter Jetv2::lower_bound_comp(Jet::SRC /**/)
{
  msg_dep_fn("lower_bound_comp");
  return Jet::lower_bound_comp(Jet::SRC::VOID);
}
Jet::Iter Jetv2::upper_bound_comp(Jet::SRC /**/)
{
  msg_dep_fn("upper_bound_comp");
  return Jet::upper_bound_comp(Jet::SRC::VOID);
}
Jet::Iter Jetv2::find(Jet::SRC /**/)
{
  msg_dep_fn("find");
  return Jet::find(Jet::SRC::VOID);
}
Jet::Iter Jetv2::end_comp()
{
  msg_dep_fn("end_comp");
  return Jet::end_comp();
}
size_t Jetv2::erase_comp(Jet::SRC /**/) {
  msg_dep_fn("erase_comp");
  return 0;
}
void Jetv2::erase_comp(Jet::Iter /*iter*/)
{
  msg_dep_fn("erase_comp");
}
void Jetv2::erase_comp(Jet::Iter /*first*/, Jet::Iter /*last*/) 
{
  msg_dep_fn("erase_comp");
}

