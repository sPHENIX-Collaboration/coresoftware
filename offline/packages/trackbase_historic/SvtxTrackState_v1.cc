#include "SvtxTrackState_v1.h"

#include <iostream>
#include <utility>   // for swap


using namespace std;

SvtxTrackState_v1::SvtxTrackState_v1(float pathlength)
  : _pathlength(pathlength)
  , _pos()
  , _mom()
  , _covar()
{
  for (int i = 0; i < 3; ++i) _pos[i] = 0.0;
  for (int i = 0; i < 3; ++i) _mom[i] = NAN;
  for (int i = 0; i < 6; ++i)
  {
    for (int j = i; j < 6; ++j)
    {
      set_error(i, j, 0.0);
    }
  }
  state_name = "UNKNOWN";
}

void SvtxTrackState_v1::identify(std::ostream &os) const
{
  os << "---SvtxTrackState_v1-------------" << endl;
  os << "pathlength: " << get_pathlength() << endl;
  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << endl;
  os << "---------------------------------" << endl;
}

float SvtxTrackState_v1::get_error(unsigned int i, unsigned int j) const
{
  return _covar[covar_index(i, j)];
}

void SvtxTrackState_v1::set_error(unsigned int i, unsigned int j, float value)
{
  _covar[covar_index(i, j)] = value;
  return;
}

unsigned int SvtxTrackState_v1::covar_index(unsigned int i, unsigned int j) const
{
  if (i > j) std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}
