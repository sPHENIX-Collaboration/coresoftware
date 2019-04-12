/**
 * @file trackbase/TrkrClusterv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrClusterv1
 */
#include "TrkrClusterv1.h"

#include <TMatrixF.h>

#include <algorithm>
#include <cmath>

TrkrClusterv1::TrkrClusterv1()
  : m_cluskey(TrkrDefs::CLUSKEYMAX)
  , m_pos()
  , m_isGlobal(true)
  , m_adc(0xFFFFFFFF)
  , m_size()
  , m_err()
{
  for (int i = 0; i < 3; ++i) m_pos[i] = NAN;

  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      setSize(i, j, NAN);
      setError(i, j, NAN);
    }
  }
}

void TrkrClusterv1::identify(std::ostream& os) const
{
  os << "---TrkrClusterv1--------------------" << std::endl;
  os << "clusid: " << getClusKey() << std::dec << std::endl;

  os << " (x,y,z) =  (" << getPosition(0);
  os << ", " << getPosition(1) << ", ";
  os << getPosition(2) << ") cm";
  if (m_isGlobal)
    os << " - global coordinates" << std::endl;
  else
    os << " - local coordinates" << std::endl;

  os << " adc = " << getAdc() << std::endl;

  os << " size phi = " << getPhiSize();
  os << " cm, size z = " << getZSize() << " cm" << std::endl;

  os << "         ( ";
  os << getSize(0, 0) << " , ";
  os << getSize(0, 1) << " , ";
  os << getSize(0, 2) << " )" << std::endl;
  os << "  size = ( ";
  os << getSize(1, 0) << " , ";
  os << getSize(1, 1) << " , ";
  os << getSize(1, 2) << " )" << std::endl;
  os << "         ( ";
  os << getSize(2, 0) << " , ";
  os << getSize(2, 1) << " , ";
  os << getSize(2, 2) << " )" << std::endl;

  os << "         ( ";
  os << getError(0, 0) << " , ";
  os << getError(0, 1) << " , ";
  os << getError(0, 2) << " )" << std::endl;
  os << "  err  = ( ";
  os << getError(1, 0) << " , ";
  os << getError(1, 1) << " , ";
  os << getError(1, 2) << " )" << std::endl;
  os << "         ( ";
  os << getError(2, 0) << " , ";
  os << getError(2, 1) << " , ";
  os << getError(2, 2) << " )" << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int TrkrClusterv1::isValid() const
{
  if (m_cluskey == TrkrDefs::CLUSKEYMAX) return 0;
  for (int i = 0; i < 3; ++i)
  {
    if (isnan(getPosition(i))) return 0;
  }
  if (m_adc == 0xFFFFFFFF) return 0;
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (isnan(getSize(i, j))) return 0;
      if (isnan(getError(i, j))) return 0;
    }
  }

  return 1;
}

void TrkrClusterv1::setSize(unsigned int i, unsigned int j, float value)
{
  m_size[covarIndex(i, j)] = value;
  return;
}

float TrkrClusterv1::getSize(unsigned int i, unsigned int j) const
{
  return m_size[covarIndex(i, j)];
}

void TrkrClusterv1::setError(unsigned int i, unsigned int j, float value)
{
  m_err[covarIndex(i, j)] = value;
  return;
}

float TrkrClusterv1::getError(unsigned int i, unsigned int j) const
{
  return m_err[covarIndex(i, j)];
}

float TrkrClusterv1::getPhiSize() const
{
  TMatrixF covar(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      covar[i][j] = getSize(i, j);
    }
  }

  float phi = -1.0 * atan2(m_pos[1], m_pos[0]);

  TMatrixF rot(3, 3);
  rot[0][0] = cos(phi);
  rot[0][1] = -sin(phi);
  rot[0][2] = 0.0;
  rot[1][0] = sin(phi);
  rot[1][1] = cos(phi);
  rot[1][2] = 0.0;
  rot[2][0] = 0.0;
  rot[2][1] = 0.0;
  rot[2][2] = 1.0;

  TMatrixF rotT(3, 3);
  rotT.Transpose(rot);

  TMatrixF trans(3, 3);
  trans = rot * covar * rotT;

  return 2.0 * sqrt(trans[1][1]);
}

float TrkrClusterv1::getZSize() const
{
  return 2.0 * sqrt(getSize(2, 2));
}

float TrkrClusterv1::getPhiError() const
{
  float rad = sqrt(m_pos[0] * m_pos[0] + m_pos[1] * m_pos[1]);
  if (rad > 0) return getRPhiError() / rad;
  return 0;
}

float TrkrClusterv1::getRPhiError() const
{
  TMatrixF covar(3, 3);
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      covar[i][j] = getError(i, j);
    }
  }

  float phi = -1.0 * atan2(m_pos[1], m_pos[0]);

  TMatrixF rot(3, 3);
  rot[0][0] = cos(phi);
  rot[0][1] = -sin(phi);
  rot[0][2] = 0.0;
  rot[1][0] = sin(phi);
  rot[1][1] = cos(phi);
  rot[1][2] = 0.0;
  rot[2][0] = 0.0;
  rot[2][1] = 0.0;
  rot[2][2] = 1.0;

  TMatrixF rotT(3, 3);
  rotT.Transpose(rot);

  TMatrixF trans(3, 3);
  trans = rot * covar * rotT;

  float rphierr = sqrt(trans[1][1]);
  if(rphierr == 0)
    {
      std::cout << "rphierr = 0 " << " x " << m_pos[0] << " y " << m_pos[1] 	<< std::endl;
      for (unsigned int i = 0; i < 3; ++i)
	{
	  for (unsigned int j = 0; j < 3; ++j)
	    {
	      std::cout << " i " << i << " j " << j << " cov " << getError(i, j) << " trans " << trans[i][j] << std::endl;
	    }
	}
    }

  return rphierr;
}

float TrkrClusterv1::getZError() const
{
  return sqrt(getError(2, 2));
}

unsigned int TrkrClusterv1::covarIndex(unsigned int i, unsigned int j) const
{
  if (i > j) std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}
