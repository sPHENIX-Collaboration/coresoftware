#include "SvtxHit_v1.h"

#include <cmath>

#include <TMatrixF.h>

using namespace std;


SvtxHit_v1::SvtxHit_v1()
  : _id(0xFFFFFFFF)
  , _layer(0xFFFFFFFF)
  , _adc(0xFFFFFFFF)
  , _e(NAN)
  , _cellid(0xFFFFFFFF)
{
}

void SvtxHit_v1::identify(ostream& os) const
{
  os << "---SvtxHit_v1--------------------" << endl;
  os << "hitid: " << get_id() << " layer: " << get_layer() << endl;
  os << "adc: " << get_adc() << " energy: " << get_e() << endl;
  os << "cellid: " << get_cellid() << endl;
  os << "---------------------------------" << endl;

  return;
}

int SvtxHit_v1::isValid() const
{
  if (_id == 0xFFFFFFFF) return 0;
  if (_layer == 0xFFFFFFFF) return 0;
  if (_adc == 0xFFFFFFFF) return 0;
  if (isnan(_e)) return 0;
  if (_cellid == 0xFFFFFFFF) return 0;
  return 1;
}
