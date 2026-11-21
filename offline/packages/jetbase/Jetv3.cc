/*!
 * \file Jetv3.cc
 * \brief Same as Jetv2.cc, plus the corrected jet pT.
 */

#include "Jetv3.h"

#include <phool/phool.h>
#include <algorithm>
#include <cmath>
#include <iostream>

class PHObject;

Jetv3::Jetv3(unsigned int n_prop)
  : _properties(n_prop, std::numeric_limits<float>::quiet_NaN())
{
}

void Jetv3::identify(std::ostream& os) const
{
  os << "---Jet v3-----------------------" << std::endl;
  os << "jetid: " << get_id() << std::endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << std::endl;
  os << " corrected pT: calib " << _pt_calib
     << ", insitu " << _pt_insitu << " (up " << _pt_insitu_up << ", down " << _pt_insitu_down << ")"
     << ", smear " << _pt_smear << " (up " << _pt_smear_up << ", down " << _pt_smear_down << ") GeV" << std::endl;

  os << " Jet Properties:";
  for (const auto& val : _properties)
  {
    os << " " << val;
  }
  os << std::endl;

  os << " Jet Components: " << _comp_ids.size() << std::endl;
  return;
}

void Jetv3::print_comp(std::ostream& os, bool single_line)
{
  for (auto iter = comp_begin(); iter != comp_end(); ++iter)
  {
    os << " (" << iter->first << "->" << static_cast<int>(iter->second) << ")";
    if (!single_line)
    {
      os << std::endl;
    }
  }
  if (single_line)
  {
    os << std::endl;
  }
}

std::vector<Jet::SRC> Jetv3::comp_src_vec()
{
  std::vector<Jet::SRC> vec{};
  auto iter = comp_begin();
  auto iter_end = comp_end();
  if (iter == iter_end)
  {
    return vec;
  }
  while (iter != iter_end)
  {
    Jet::SRC src = iter->first;
    vec.push_back(src);
    iter = comp_end(src);
  }
  return vec;
}

std::map<Jet::SRC, size_t> Jetv3::comp_src_sizemap()
{
  std::map<Jet::SRC, size_t> sizemap{};
  auto iter = comp_begin();
  auto iter_end = comp_end();
  if (iter == iter_end)
  {
    return sizemap;
  }
  while (iter != iter_end)
  {
    Jet::SRC src = iter->first;
    auto iter_ub = comp_end(src);
    sizemap.insert(std::make_pair(src, static_cast<size_t>(iter_ub - iter)));
    iter = iter_ub;
  }
  return sizemap;
}

void Jetv3::Reset()
{
  _id = 0xFFFFFFFF;
  std::fill(std::begin(_mom), std::end(_mom), std::numeric_limits<float>::quiet_NaN());
  _e = std::numeric_limits<float>::quiet_NaN();
  _emcal_frac = std::numeric_limits<float>::quiet_NaN();
  _ihcal_frac = std::numeric_limits<float>::quiet_NaN();
  _ohcal_frac = std::numeric_limits<float>::quiet_NaN();
  _pt_calib = std::numeric_limits<float>::quiet_NaN();
  _pt_insitu = std::numeric_limits<float>::quiet_NaN();
  _pt_insitu_up = std::numeric_limits<float>::quiet_NaN();
  _pt_insitu_down = std::numeric_limits<float>::quiet_NaN();
  _pt_smear = std::numeric_limits<float>::quiet_NaN();
  _pt_smear_up = std::numeric_limits<float>::quiet_NaN();
  _pt_smear_down = std::numeric_limits<float>::quiet_NaN();
  // also reset (unlike Jetv2): a jet slot reused by the TClonesArray must not carry these over
  _isCalib = 0;
  _is_sorted = false;
  _comp_ids.clear();
  _properties.clear();
}

int Jetv3::isValid() const
{
  if (_id == 0xFFFFFFFF)
  {
    return 0;
  }
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
  if (_comp_ids.empty())
  {
    return 0;
  }
  return 1;
}

PHObject* Jetv3::CloneMe() const
{
  Jet* jet = new Jetv3(*this);
  return jet;
}

float Jetv3::get_p() const
{
  return std::sqrt((get_px() * get_px()) + (get_py() * get_py()) + (get_pz() * get_pz()));
}

float Jetv3::get_pt() const
{
  return std::sqrt((get_px() * get_px()) + (get_py() * get_py()));
}

float Jetv3::get_et() const
{
  return get_pt() / get_p() * get_e();
}

float Jetv3::get_eta() const
{
  return std::asinh(get_pz() / get_pt());
}

float Jetv3::get_phi() const
{
  return std::atan2(get_py(), get_px());
}

float Jetv3::get_mass() const
{
  // follow CLHEP convention and return negative mass if E^2 - p^2 < 0
  float mass2 = get_mass2();
  if (mass2 < 0)
  {
    return -1 * std::sqrt(std::fabs(mass2));
  }
  return std::sqrt(mass2);
}

float Jetv3::get_mass2() const
{
  float p2 = (get_px() * get_px()) + (get_py() * get_py()) + (get_pz() * get_pz());
  return (get_e() * get_e()) - p2;
}

size_t Jetv3::num_comp(Jet::SRC iSRC)
{
  if (iSRC == Jet::SRC::VOID)
  {
    return (comp_end() - comp_begin());
  }
  return (comp_end(iSRC) - comp_begin(iSRC));
}

void Jetv3::insert_comp(SRC iSRC, unsigned int compid)
{
  _is_sorted = false;
  _comp_ids.emplace_back(iSRC, compid);
}
void Jetv3::insert_comp(SRC iSRC, unsigned int compid, bool /*skip_flag*/)
{
  _comp_ids.emplace_back(iSRC, compid);
}

// same as above, but add whole vector
void Jetv3::insert_comp(Jet::TYPE_comp_vec& added)
{
  _is_sorted = false;
  _comp_ids.insert(_comp_ids.end(), added.begin(), added.end());
}
void Jetv3::insert_comp(Jet::TYPE_comp_vec& added, bool /*skip_nosort_flag_set*/)
{
  _comp_ids.insert(_comp_ids.end(), added.begin(), added.end());
}

void Jetv3::ensure_sorted()
{
  if (!_is_sorted)
  {
    std::sort(_comp_ids.begin(), _comp_ids.end(), CompareSRC());
    _is_sorted = true;
  }
}

Jetv3::ITER_comp_vec Jetv3::comp_begin(Jet::SRC iSRC)
{
  ensure_sorted();
  return std::lower_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

Jetv3::ITER_comp_vec Jetv3::comp_end(Jet::SRC iSRC)
{
  ensure_sorted();
  return std::upper_bound(_comp_ids.begin(), _comp_ids.end(), iSRC, CompareSRC());
}

void Jetv3::msg_dep_fn(const std::string& fn_name)
{
  std::cout << " warning: Method Jet::" << fn_name << "() deprecated in Jetv3" << std::endl;
}

bool Jetv3::has_property(Jet::PROPERTY /*prop_id*/) const
{
  msg_dep_fn("has_property");
  return false;
}

void Jetv3::print_property(std::ostream& /**/) const
{
  msg_dep_fn("print_property");
}

bool Jetv3::empty_comp() const
{
  msg_dep_fn("empty_comp");
  return true;
}

size_t Jetv3::count_comp(Jet::SRC /**/) const
{
  msg_dep_fn("count_comp");
  return 0;
}

Jet::ConstIter Jetv3::begin_comp() const
{
  msg_dep_fn("begin_comp");
  return Jet::begin_comp();
}
Jet::ConstIter Jetv3::lower_bound_comp(Jet::SRC /**/) const
{
  msg_dep_fn("lower_bound_comp");
  return Jet::lower_bound_comp(Jet::SRC::VOID);
}
Jet::ConstIter Jetv3::upper_bound_comp(Jet::SRC /**/) const
{
  msg_dep_fn("upper_bound_comp");
  return Jet::upper_bound_comp(Jet::SRC::VOID);
}
Jet::ConstIter Jetv3::find(Jet::SRC /**/) const
{
  msg_dep_fn("find");
  return Jet::find(Jet::SRC::VOID);
}
Jet::ConstIter Jetv3::end_comp() const
{
  msg_dep_fn("end_comp");
  return Jet::end_comp();
}
Jet::Iter Jetv3::begin_comp()
{
  msg_dep_fn("begin_comp");
  return Jet::begin_comp();
}
Jet::Iter Jetv3::lower_bound_comp(Jet::SRC /**/)
{
  msg_dep_fn("lower_bound_comp");
  return Jet::lower_bound_comp(Jet::SRC::VOID);
}
Jet::Iter Jetv3::upper_bound_comp(Jet::SRC /**/)
{
  msg_dep_fn("upper_bound_comp");
  return Jet::upper_bound_comp(Jet::SRC::VOID);
}
Jet::Iter Jetv3::find(Jet::SRC /**/)
{
  msg_dep_fn("find");
  return Jet::find(Jet::SRC::VOID);
}
Jet::Iter Jetv3::end_comp()
{
  msg_dep_fn("end_comp");
  return Jet::end_comp();
}
size_t Jetv3::erase_comp(Jet::SRC /**/)
{
  msg_dep_fn("erase_comp");
  return 0;
}
void Jetv3::erase_comp(Jet::Iter /*iter*/)
{
  msg_dep_fn("erase_comp");
}
void Jetv3::erase_comp(Jet::Iter /*first*/, Jet::Iter /*last*/)
{
  msg_dep_fn("erase_comp");
}
