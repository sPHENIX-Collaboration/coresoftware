#include "TrkrClusterv1.h"

#include <TMatrixF.h>

#include <algorithm>
#include <cmath>

using namespace std;

ClassImp(TrkrClusterv1);

TrkrClusterv1::TrkrClusterv1()
  : cluskey_(TrkrDefs::CLUSKEYMAX)
  , pos_()
  , is_global_(true)
  , adc_(0xFFFFFFFF)
  , size_()
  , err_()
{
  for (int i = 0; i < 3; ++i) pos_[i] = NAN;

  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      SetSize(i, j, NAN);
      SetError(i, j, NAN);
    }
  }
}

void TrkrClusterv1::identify(ostream& os) const
{
  os << "---TrkrClusterv1--------------------" << endl;
  os << "clusid: 0x" << std::hex << GetClusKey() << std::dec << endl;

  os << " (x,y,z) =  (" << GetPosition(0);
  os << ", " << GetPosition(1) << ", ";
  os << GetPosition(2) << ") cm";
  if (is_global_)
    os << " - global coordinates" << endl;
  else
    os << " - local coordinates" << endl;

  os << " adc = " << GetAdc() << endl;

  os << " size phi = " << GetPhiSize();
  os << " cm, size z = " << GetZSize() << " cm" << endl;

  os << "         ( ";
  os << GetSize(0, 0) << " , ";
  os << GetSize(0, 1) << " , ";
  os << GetSize(0, 2) << " )" << endl;
  os << "  size = ( ";
  os << GetSize(1, 0) << " , ";
  os << GetSize(1, 1) << " , ";
  os << GetSize(1, 2) << " )" << endl;
  os << "         ( ";
  os << GetSize(2, 0) << " , ";
  os << GetSize(2, 1) << " , ";
  os << GetSize(2, 2) << " )" << endl;

  os << "         ( ";
  os << GetError(0, 0) << " , ";
  os << GetError(0, 1) << " , ";
  os << GetError(0, 2) << " )" << endl;
  os << "  err  = ( ";
  os << GetError(1, 0) << " , ";
  os << GetError(1, 1) << " , ";
  os << GetError(1, 2) << " )" << endl;
  os << "         ( ";
  os << GetError(2, 0) << " , ";
  os << GetError(2, 1) << " , ";
  os << GetError(2, 2) << " )" << endl;

  os << endl;
  os << "-----------------------------------------------" << endl;

  return;
}

int TrkrClusterv1::isValid() const
{
  if (cluskey_ == TrkrDefs::CLUSKEYMAX) return 0;
  for (int i = 0; i < 3; ++i)
  {
    if (isnan(GetPosition(i))) return 0;
  }
  if (adc_ == 0xFFFFFFFF) return 0;
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (isnan(GetSize(i, j))) return 0;
      if (isnan(GetError(i, j))) return 0;
    }
  }

  return 1;
}

void TrkrClusterv1::SetSize(unsigned int i, unsigned int j, float value)
{
  size_[CovarIndex(i, j)] = value;
  return;
}

float TrkrClusterv1::GetSize(unsigned int i, unsigned int j) const
{
  return size_[CovarIndex(i, j)];
}

void TrkrClusterv1::SetError(unsigned int i, unsigned int j, float value)
{
  err_[CovarIndex(i, j)] = value;
  return;
}

float TrkrClusterv1::GetError(unsigned int i, unsigned int j) const
{
  return err_[CovarIndex(i, j)];
}

float TrkrClusterv1::GetPhiSize() const
{
  TMatrixF COVAR(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      COVAR[i][j] = GetSize(i, j);
    }
  }

  float phi = -1.0 * atan2(pos_[1], pos_[0]);

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

float TrkrClusterv1::GetZSize() const
{
  return 2.0 * sqrt(GetSize(2, 2));
}

float TrkrClusterv1::GetPhiError() const
{
  float rad = sqrt(pos_[0] * pos_[0] + pos_[1] * pos_[1]);
  if (rad > 0) return GetRPhiError() / rad;
  return 0;
}

float TrkrClusterv1::GetRPhiError() const
{
  TMatrixF COVAR(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      COVAR[i][j] = GetError(i, j);
    }
  }

  float phi = -1.0 * atan2(pos_[1], pos_[0]);

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

float TrkrClusterv1::GetZError() const
{
  return sqrt(GetError(2, 2));
}

unsigned int TrkrClusterv1::CovarIndex(unsigned int i, unsigned int j) const
{
  if (i > j) std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}
