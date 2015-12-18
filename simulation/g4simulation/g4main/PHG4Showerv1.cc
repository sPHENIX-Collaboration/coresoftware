#include "PHG4Showerv1.h"

#include "PHG4HitDefs.h"

#include <cmath>
#include <iostream>

using namespace std;

ClassImp(PHG4Showerv1);

PHG4Showerv1::PHG4Showerv1()
  : _id(0xFFFFFFFF), _primary_id(-1), _parent_shower_id(0),
    _pos(), _covar(), _edep(), _eion(),
    _light_yield(), _g4particle_ids(), _g4hit_ids() {

  for (int i = 0; i < 3; ++i)
    _pos[i] = NAN;

  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_covar(i, j, NAN);
    }
  }
}

void PHG4Showerv1::identify(ostream &os) const {
  os << "---PHG4Showerv1-------------------------------" << endl;
  os << "id: " << get_id() << endl;
  os << "primary_id: " << get_primary_id() << endl;
  os << "x: " << get_x() << endl;
  os << "y: " << get_y() << endl;
  os << "z: " << get_z() << endl;

  os << "         ( ";
  os << get_covar(0, 0) << " , ";
  os << get_covar(0, 1) << " , ";
  os << get_covar(0, 2) << " )" << endl;
  os << " covar = ( ";
  os << get_covar(1, 0) << " , ";
  os << get_covar(1, 1) << " , ";
  os << get_covar(1, 2) << " )" << endl;
  os << "         ( ";
  os << get_covar(2, 0) << " , ";
  os << get_covar(2, 1) << " , ";
  os << get_covar(2, 2) << " )" << endl;

  os << "VOLUME ID : edep eion light_yield" << endl;
  for (std::map<int, float>::const_iterator iter = _edep.begin();
       iter != _edep.end(); ++iter) {
    int volid = iter->first;
    os << volid << " : " << get_edep(volid) << " " << get_eion(volid) << " "
       << get_light_yield(volid) << endl;
  }

  os << "G4Particle IDs" << endl;
  for (std::set<int>::const_iterator iter = _g4particle_ids.begin();
       iter != _g4particle_ids.end(); ++iter) {
    os << *iter << " ";
  }
  os << endl;

  os << "G4Hit IDs" << endl;
  for (std::map<int,std::set<PHG4HitDefs::keytype> >::const_iterator iter = _g4hit_ids.begin();
       iter != _g4hit_ids.end();
       ++iter) {
    for (std::set<PHG4HitDefs::keytype>::const_iterator jter = iter->second.begin();
	 jter != iter->second.end(); ++jter) {
      os << *jter << " ";
    }
  }
  os << endl;

  os << "-----------------------------------------------" << endl;

  return;
}

int PHG4Showerv1::isValid() const {
  if (_id == 0xFFFFFFFF)
    return 0;
  if (_primary_id == -1)
    return 0;
  for (int i = 0; i < 3; ++i) {
    if (isnan(_pos[i]))
      return 0;
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      if (isnan(get_covar(i, j)))
        return 0;
    }
  }
  return 1;
}

void PHG4Showerv1::set_covar(unsigned int i, unsigned int j, float value) {
  _covar[covar_index(i, j)] = value;
  return;
}

float PHG4Showerv1::get_covar(unsigned int i, unsigned int j) const {
  return _covar[covar_index(i, j)];
}

unsigned int PHG4Showerv1::covar_index(unsigned int i, unsigned int j) const {
  if (i > j)
    std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

float PHG4Showerv1::get_edep(int volume) const {
  std::map<int, float>::const_iterator citer =
      _edep.find(volume);
  if (citer == _edep.end())
    return 0.0;
  return citer->second;
}

float PHG4Showerv1::get_eion(int volume) const {
  std::map<int, float>::const_iterator citer =
      _eion.find(volume);
  if (citer == _eion.end())
    return 0.0;
  return citer->second;
}

float PHG4Showerv1::get_light_yield(int volume) const {
  std::map<int, float>::const_iterator citer =
      _light_yield.find(volume);
  if (citer == _light_yield.end())
    return 0.0;
  return citer->second;
}
