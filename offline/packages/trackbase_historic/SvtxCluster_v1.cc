#include "SvtxCluster_v1.h"

#include <TMatrixF.h>

#include <algorithm>
#include <cmath>

using namespace std;

ClassImp(SvtxCluster_v1);

SvtxCluster_v1::SvtxCluster_v1()
  : _id(0xFFFFFFFF)
  , _layer(0xFFFFFFFF)
  , _pos()
  , _e(NAN)
  , _adc(0xFFFFFFFF)
  , _size()
  , _err()
  , _hit_ids()
{
  for (int i = 0; i < 3; ++i) _pos[i] = NAN;

  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      set_size(i, j, NAN);
      set_error(i, j, NAN);
    }
  }
}

void SvtxCluster_v1::identify(ostream& os) const
{
  os << "---SvtxCluster_v1--------------------" << endl;
  os << "clusid: " << get_id() << " layer: " << get_layer() << endl;

  os << " (x,y,z) =  (" << get_position(0);
  os << ", " << get_position(1) << ", ";
  os << get_position(2) << ") cm" << endl;

  os << " e = " << get_e() << " adc = " << get_adc() << endl;

  os << " size phi = " << get_phi_size();
  os << " cm, size z = " << get_z_size() << " cm" << endl;

  os << "         ( ";
  os << get_size(0, 0) << " , ";
  os << get_size(0, 1) << " , ";
  os << get_size(0, 2) << " )" << endl;
  os << "  size = ( ";
  os << get_size(1, 0) << " , ";
  os << get_size(1, 1) << " , ";
  os << get_size(1, 2) << " )" << endl;
  os << "         ( ";
  os << get_size(2, 0) << " , ";
  os << get_size(2, 1) << " , ";
  os << get_size(2, 2) << " )" << endl;

  os << "         ( ";
  os << get_error(0, 0) << " , ";
  os << get_error(0, 1) << " , ";
  os << get_error(0, 2) << " )" << endl;
  os << "  err  = ( ";
  os << get_error(1, 0) << " , ";
  os << get_error(1, 1) << " , ";
  os << get_error(1, 2) << " )" << endl;
  os << "         ( ";
  os << get_error(2, 0) << " , ";
  os << get_error(2, 1) << " , ";
  os << get_error(2, 2) << " )" << endl;

  os << " list of hits ids: ";
  for (ConstHitIter iter = begin_hits(); iter != end_hits(); ++iter)
  {
    os << *iter << " ";
  }
  os << endl;
  os << "-----------------------------------------------" << endl;

  return;
}

int SvtxCluster_v1::isValid() const
{
  if (_id == 0xFFFFFFFF) return 0;
  if (_layer == 0xFFFFFFFF) return 0;
  for (int i = 0; i < 3; ++i)
  {
    if (isnan(_pos[i])) return 0;
  }
  if (isnan(_e)) return 0;
  if (_adc == 0xFFFFFFFF) return 0;
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (isnan(get_size(i, j))) return 0;
      if (isnan(get_error(i, j))) return 0;
    }
  }
  if (_hit_ids.empty()) return 0;

  return 1;
}

void SvtxCluster_v1::set_size(unsigned int i, unsigned int j, float value)
{
  _size[covar_index(i, j)] = value;
  return;
}

float SvtxCluster_v1::get_size(unsigned int i, unsigned int j) const
{
  return _size[covar_index(i, j)];
}

void SvtxCluster_v1::set_error(unsigned int i, unsigned int j, float value)
{
  _err[covar_index(i, j)] = value;
  return;
}

float SvtxCluster_v1::get_error(unsigned int i, unsigned int j) const
{
  return _err[covar_index(i, j)];
}

float SvtxCluster_v1::get_phi_size() const
{
  TMatrixF COVAR(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      COVAR[i][j] = get_size(i, j);
    }
  }

  float phi = -1.0 * atan2(_pos[1], _pos[0]);

  TMatrixF ROT(3, 3);
  ROT[0][0] = cos(phi);
  ROT[0][1] = -sin(phi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(phi);
  ROT[1][1] = cos(phi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0;
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;

  TMatrixF ROT_T(3, 3);
  ROT_T.Transpose(ROT);

  TMatrixF TRANS(3, 3);
  TRANS = ROT * COVAR * ROT_T;

  return 2.0 * sqrt(TRANS[1][1]);
}

float SvtxCluster_v1::get_z_size() const
{
  return 2.0 * sqrt(get_size(2, 2));
}

float SvtxCluster_v1::get_phi_error() const
{
  float rad = sqrt(_pos[0] * _pos[0] + _pos[1] * _pos[1]);
  if (rad > 0) return get_rphi_error() / rad;
  return 0;
}

float SvtxCluster_v1::get_rphi_error() const
{
  TMatrixF COVAR(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      COVAR[i][j] = get_error(i, j);
    }
  }

  float phi = -1.0 * atan2(_pos[1], _pos[0]);

  TMatrixF ROT(3, 3);
  ROT[0][0] = cos(phi);
  ROT[0][1] = -sin(phi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(phi);
  ROT[1][1] = cos(phi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0;
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;

  TMatrixF ROT_T(3, 3);
  ROT_T.Transpose(ROT);

  TMatrixF TRANS(3, 3);
  TRANS = ROT * COVAR * ROT_T;

  return sqrt(TRANS[1][1]);
}

float SvtxCluster_v1::get_z_error() const
{
  return sqrt(get_error(2, 2));
}

unsigned int SvtxCluster_v1::covar_index(unsigned int i, unsigned int j) const
{
  if (i > j) std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}
