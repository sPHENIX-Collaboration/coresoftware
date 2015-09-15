#include "SvtxHit.h"

#include <cmath>

#include <TMatrixF.h>

using namespace std;

ClassImp(SvtxHit);

SvtxHit::SvtxHit()
  : _id(0xFFFFFFFF),
    _layer(0xFFFFFFFF),
    _adc(0xFFFFFFFF),
    _e(NAN),
    _cellid(0xFFFFFFFF)
{}

void SvtxHit::identify(ostream& os) const {
  os << "---SvtxHit-----------------------" << endl;
  os << "hitid: " << get_id() << " layer: "<< get_layer() << endl;
  os << "adc: " << get_adc() << " energy: "<< get_e() << endl;
  os << "cellid: " << get_cellid() << endl;
  os << "---------------------------------" << endl;
  
  return;  
}

void SvtxHit::Reset() {
  *this = SvtxHit();
}

int SvtxHit::IsValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  if (_layer == 0xFFFFFFFF) return 0;
  if (_adc == 0xFFFFFFFF) return 0;
  if (isnan(_e)) return 0;
  if (_cellid == 0xFFFFFFFF) return 0;
  return 1;
}
