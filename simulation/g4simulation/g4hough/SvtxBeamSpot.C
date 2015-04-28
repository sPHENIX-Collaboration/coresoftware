#include "SvtxBeamSpot.h"

#include <cmath>

using namespace std;

ClassImp(SvtxBeamSpot);

SvtxBeamSpot::SvtxBeamSpot()
  : _pos(),
    _err() {
  
  for (int i = 0; i < 2; ++i) _pos[i] = NAN;
  for (int i = 0; i < 2; ++i) _err[i] = new float[i+1];

  for (int j = 0; j < 2; ++j) {
    for (int i = j; i < 2; ++i) {
      set_error(i,j,NAN);
    }
  } 
}

SvtxBeamSpot::SvtxBeamSpot(const SvtxBeamSpot &beamspot) :
  _pos(),
  _err() {
  
  for (int i=0; i<2; ++i) _pos[i] = beamspot.get_position(i);    
  for (int i=0; i<2; ++i) _err[i] = new float[i+1];
  for (int j=0; j<2; ++j) {
    for (int i=j; i<2; ++i) {
      set_error(i,j,beamspot.get_error(i,j));
    }
  } 
}

SvtxBeamSpot& SvtxBeamSpot::operator=(const SvtxBeamSpot &beamspot) {
  Reset();
  for (int i = 0; i < 2; ++i) _pos[i] = beamspot.get_position(i);
  for (int j = 0; j < 2; ++j) {
    for (int i = j; i < 2; ++i) {
      set_error(i,j,beamspot.get_error(i,j));
    }
  } 

  return *this;
}

SvtxBeamSpot::~SvtxBeamSpot(){
  for (int i=0; i<2; ++i) delete[] _err[i];
}

void SvtxBeamSpot::identify(ostream& os) const {
  os << "---SvtxBeamSpot--------------------------------" << endl;
  
  os << " (x,y) =  (" << get_position(0);
  os << ", " << get_position(1) << ") cm" << endl;

  os << "  err = ( ";
  os << get_error(0,0) << " , ";
  os << get_error(0,1) << " )" << endl; 
  os << "        ( ";
  os << get_error(1,0) << " , ";
  os << get_error(1,1) << " )" << endl;
  os << "-----------------------------------------------" << endl;
  
  return;  
}

void SvtxBeamSpot::Reset() {
  for (int i = 0; i < 2; ++i) _pos[i] = NAN;
  for (int j = 0; j < 2; ++j) {
    for (int i = j; i < 2; ++i) {
      set_error(i,j,NAN);
    }
  } 
}

int SvtxBeamSpot::IsValid() const {
  for (int i = 0; i < 2; ++i) {
    if (isnan(_pos[i])) return 0;
  }
  for (int j = 0; j < 2; ++j) {
    for (int i = j; i < 2; ++i) {
      if (isnan(get_error(i,j))) return 0;
    }
  } 
  
  return 1;
}

void SvtxBeamSpot::set_error(int i, int j, float value) {
  if (j > i) set_error(j,i,value);
  else _err[i][j] = value;
  return;
}

float SvtxBeamSpot::get_error(int i, int j) const {
  if (j > i) return get_error(j,i);
  return _err[i][j];
}
