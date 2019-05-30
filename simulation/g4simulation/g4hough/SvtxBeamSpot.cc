#include "SvtxBeamSpot.h"

#include <cmath>
#include <utility>  // for swap

using namespace std;

SvtxBeamSpot::SvtxBeamSpot()
  : _pos(),
    _err() { 
  for (int i = 0; i < 2; ++i) _pos[i] = NAN;
  for (int j = 0; j < 2; ++j) {
    for (int i = j; i < 2; ++i) {
      set_error(i,j,NAN);
    }
  } 
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

int SvtxBeamSpot::isValid() const {
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
  _err[covar_index(i,j)] = value;
  return;
}

float SvtxBeamSpot::get_error(int i, int j) const {
  return _err[covar_index(i,j)];
}

unsigned int SvtxBeamSpot::covar_index(unsigned int i, unsigned int j) const {
  if (i>j) std::swap(i,j);
  return i+1+(j+1)*(j)/2-1;
}

