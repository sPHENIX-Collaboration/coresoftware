#include "PHG4Shower_v1.h"

#include <cmath>
#include <iostream>

using namespace std;

ClassImp(PHG4Shower_v1);

PHG4Shower_v1::PHG4Shower_v1()
  : _id(0xFFFFFFFF),
    _primary_id(-1),
    _pos(),
    _covar(),
    _edep(),
    _eion(),
    _light_yield() {

  for (int i = 0; i < 3; ++i) _pos[i] = NAN;  

  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_covar(i,j,NAN);
    }
  }   
}

PHG4Shower_v1::~PHG4Shower_v1(){}

void PHG4Shower_v1::identify(ostream& os) const {
  os << "---PHG4Shower_v1-------------------------------" << endl;
  os << "id: " << get_id() << endl;
  os << "primary_id: " << get_primary_id() << endl;
  os << "x: " << get_x() << endl;
  os << "y: " << get_y() << endl;
  os << "z: " << get_z() << endl;

  os << "         ( ";
  os << get_covar(0,0) << " , ";
  os << get_covar(0,1) << " , ";
  os << get_covar(0,2) << " )" << endl; 
  os << " covar = ( ";
  os << get_covar(1,0) << " , ";
  os << get_covar(1,1) << " , ";
  os << get_covar(1,2) << " )" << endl;
  os << "         ( ";
  os << get_covar(2,0) << " , ";
  os << get_covar(2,1) << " , ";
  os << get_covar(2,2) << " )" << endl;

  os << "VOLUME ID : edep eion light_yield" << endl;
  for (std::map<PHG4Shower::VOLUME,float>::const_iterator iter = _edep.begin();
       iter != _edep.end();
       ++iter) {
    PHG4Shower::VOLUME volid = iter->first;
    cout << volid << " : "
	 << get_edep(volid) << " "
      	 << get_eion(volid) << " "
      	 << get_light_yield(volid) << endl;
  }
  
  os << "-----------------------------------------------" << endl;
  
  return;  
}

PHG4Shower* PHG4Shower_v1::Clone() const {
  PHG4Shower_v1* shower = new PHG4Shower_v1(*this);
  return shower;
}

void PHG4Shower_v1::Reset() {
  _id = 0xFFFFFFFF;
  _primary_id = -1;
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;  
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      set_covar(i,j,NAN);
    }
  }
  _edep.clear();
  _eion.clear();
  _light_yield.clear();
}

int PHG4Shower_v1::isValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  if (_primary_id == -1) return 0;
  for (int i = 0; i < 3; ++i) {
    if (isnan(_pos[i])) return 0;
  }
  for (int j = 0; j < 3; ++j) {
    for (int i = j; i < 3; ++i) {
      if (isnan(get_covar(i,j))) return 0;
    }
  }
  if (_edep.empty()) return 0;
  if (_eion.empty()) return 0;
  if (_light_yield.empty()) return 0;
  return 1;
}

void PHG4Shower_v1::set_covar(unsigned int i, unsigned int j, float value) {
  _covar[covar_index(i,j)] = value;
  return;
}

float PHG4Shower_v1::get_covar(unsigned int i, unsigned int j) const {
  return _covar[covar_index(i,j)];
}

unsigned int PHG4Shower_v1::covar_index(unsigned int i, unsigned int j) const {
  if (i>j) std::swap(i,j);
  return i+1+(j+1)*(j)/2-1;
}

float PHG4Shower_v1::get_edep(PHG4Shower::VOLUME volume) const {
  std::map<PHG4Shower::VOLUME,float>::const_iterator citer = _edep.find(volume);
  if (citer == _edep.end()) return NAN;
  return citer->second;
}

float PHG4Shower_v1::get_eion(PHG4Shower::VOLUME volume) const {
  std::map<PHG4Shower::VOLUME,float>::const_iterator citer = _eion.find(volume);
  if (citer == _eion.end()) return NAN;
  return citer->second;
}

float PHG4Shower_v1::get_light_yield(PHG4Shower::VOLUME volume) const {
  std::map<PHG4Shower::VOLUME,float>::const_iterator citer = _light_yield.find(volume);
  if (citer == _light_yield.end()) return NAN;
  return citer->second;
}
