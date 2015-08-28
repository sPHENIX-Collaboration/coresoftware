#include "BbcVertex.h"

#include <cmath>

using namespace std;

ClassImp(BbcVertex);

BbcVertex::BbcVertex()
  : _id(0xFFFFFFFF),
    _t0(NAN),
    _t0_err(NAN),
    _z(NAN),
    _z_err(NAN) {
}

BbcVertex::~BbcVertex(){}

void BbcVertex::identify(ostream& os) const {
  os << "---BbcVertex-----------------------------------" << endl;
  os << "vertexid: " << get_id() << endl;
  os << " t0 = " << get_t0() << " +/- " << get_t0_err() << endl;
  os << " z =  " << get_z() << " +/- " << get_z_err() << endl;
  os << "-----------------------------------------------" << endl;
  
  return;  
}

void BbcVertex::Reset() {
  _id = 0xFFFFFFFF;
  _t0 = NAN;
  _t0_err = NAN;
  _z = NAN;
  _z_err = NAN;
}

int BbcVertex::IsValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  if (isnan(_t0)) return 0;
  if (isnan(_t0_err)) return 0;
  if (isnan(_z)) return 0;
  if (isnan(_z_err)) return 0;
  
  return 1;
}
