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

SvtxHit::SvtxHit(const SvtxHit &hit) :
  _id(hit.get_id()),
  _layer(hit.get_layer()),
  _adc(hit.get_adc()),
  _e(hit.get_e()),
  _cellid(hit.get_cellid())
{}

SvtxHit& SvtxHit::operator=(const SvtxHit &hit) {
  Reset();
  _id = hit.get_id();
  _layer = hit.get_layer();
  _adc = hit.get_adc();
  _e = hit.get_e();
  _cellid = hit.get_cellid();

  return *this;
}

SvtxHit::~SvtxHit(){}

void SvtxHit::identify(ostream& os) const {
  os << "---SvtxHit-----------------------" << endl;
  os << "hitid: " << get_id() << " layer: "<< get_layer() << endl;
  os << "adc: " << get_adc() << " energy: "<< get_e() << endl;
  os << "cellid: " << get_cellid() << endl;
  os << "---------------------------------" << endl;
  
  return;  
}

void SvtxHit::Reset() {
  _id     = 0xFFFFFFFF;
  _layer  = 0xFFFFFFFF;
  _adc    = 0xFFFFFFFF;
  _e      = NAN;
  _cellid = 0xFFFFFFFF;
}

int SvtxHit::IsValid() const {
  if (_id == 0xFFFFFFFF) return 0;
  if (_layer == 0xFFFFFFFF) return 0;
  if (_adc == 0xFFFFFFFF) return 0;
  if (isnan(_e)) return 0;
  if (_cellid == 0xFFFFFFFF) return 0;
  return 1;
}
