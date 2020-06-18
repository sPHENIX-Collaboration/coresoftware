#include "SvtxTrackState_v1.h"

#include <iostream>
#include <utility>   // for swap


using namespace std;

namespace
{

  // square convenience function
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  // get unique index in cov. matrix array from i and j
  inline unsigned int covar_index(unsigned int i, unsigned int j)
  {
    if (i > j) std::swap(i, j);
    return i + 1 + (j + 1) * (j) / 2 - 1;
  }

}

SvtxTrackState_v1::SvtxTrackState_v1(float pathlength)
  : _pathlength(pathlength)
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

float SvtxTrackState_v1::get_phi_error() const
{
  const float r = std::sqrt( square(_pos[0]) + square(_pos[1]));
  if (r > 0) return get_rphi_error() / r;
  return 0;
}

float SvtxTrackState_v1::get_rphi_error() const
{
  const auto phi = -std::atan2(_pos[1], _pos[0] );
  const auto cosphi = std::cos(phi);
  const auto sinphi = std::sin(phi);
  return std::sqrt(
    square(sinphi)*get_error(0,0) +
    square(cosphi)*get_error(1,1) +
    2.*cosphi*sinphi*get_error(0,1));
}

float SvtxTrackState_v1::get_z_error() const
{ return std::sqrt(get_error(2, 2)); }
