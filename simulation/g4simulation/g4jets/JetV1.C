/*!
 * \file JetV1.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "JetV1.h"

#include <cmath>

using namespace std;

ClassImp(JetV1);

JetV1::JetV1()
  : _id(0xFFFFFFFF),
    _mom(),
    _e(NAN),
    _comp_ids(),
    _property_map() {
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
}

void JetV1::identify(ostream& os) const {
  os << "---Jet V1-----------------------" << endl;
  os << "jetid: " << get_id() << endl;
  os << " (px,py,pz,e) =  (" << get_px() << ", " << get_py() << ", ";
  os << get_pz() << ", " << get_e() << ") GeV" << endl;
  print_property(os);
  for (ConstIter citer = begin_comp(); citer != end_comp(); ++citer) {
    cout << citer->first << " -> " << citer->second << endl;
  }
  os << "-----------------------------------------------" << endl;

  return;
}

void JetV1::Reset() {
  _id = 0xFFFFFFFF;
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
  _e = NAN;
  _comp_ids.clear();
  _property_map.clear();
}

int JetV1::isValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  for (int i = 0; i < 3; ++i) {
    if (isnan(_mom[i])) return 0;
  }
  if (isnan(_e)) return 0;
  if (_comp_ids.empty()) return 0;
  return 1;
}

float JetV1::get_p() const {
  return sqrt(get_px()*get_px()+get_py()*get_py()+get_pz()*get_pz());
}

float JetV1::get_pt() const {
  return sqrt(get_px()*get_px()+get_py()*get_py());
}

float JetV1::get_et() const {
  return get_pt()/get_p()*get_e();
}

float JetV1::get_eta() const {
  return asinh(get_pz()/get_pt());
}

float JetV1::get_phi() const {
  return atan2(get_py(),get_px());
}

float JetV1::get_mass() const {
  float p2 = get_px()*get_px()+get_py()*get_py()+get_pz()*get_pz();
  return sqrt(get_e()*get_e()-p2);
}

bool JetV1::has_property(Jet::PROPERTY prop_id) const {
  typ_property_map::const_iterator citer = _property_map.find(prop_id); 
  if (citer==_property_map.end()) return false;
  else return true;
}

float JetV1::get_property(Jet::PROPERTY prop_id) const {
  typ_property_map::const_iterator citer = _property_map.find(prop_id);
  if (citer==_property_map.end()) return NAN;
  else return citer->second;
}

void JetV1::set_property(Jet::PROPERTY prop_id, float value) {
  _property_map[prop_id] = value;
}

void JetV1::print_property(ostream& os) const {
  for (typ_property_map::const_iterator citer = _property_map.begin();
      citer != _property_map.end(); ++citer) {

    cout << " "; //indent

    switch (citer->first) {
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
