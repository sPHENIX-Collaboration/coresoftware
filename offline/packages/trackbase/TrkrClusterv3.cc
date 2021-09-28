/**
 * @file trackbase/TrkrClusterv3.cc
 * @author J. Osborn
 * @date October 2021
 * @brief Implementation of TrkrClusterv3
 */
#include "TrkrClusterv3.h"

#include <cmath>
#include <utility>          // for swap

namespace
{
  // square convenience function
  template<class T> inline constexpr T square( const T& x ) 
    { return x*x; }

  // get unique index in cov. matrix array from i and j
  inline unsigned int covarIndex(unsigned int i, unsigned int j)
  {
    if (i > j) std::swap(i, j);
    return i + 1 + (j + 1) * (j) / 2 - 1;
  }

 // rotate size or covariant matrix to polar coordinates and return the phi component
 template<float (TrkrClusterv3::*accessor)(unsigned int, unsigned int) const>
    float rotate( const TrkrClusterv3* cluster )
  {
    const auto phi = -std::atan2(cluster->getY(), cluster->getX());
    const auto cosphi = std::cos(phi);
    const auto sinphi = std::sin(phi);

    return
      square(sinphi)*(cluster->*accessor)(0,0) +
      square(cosphi)*(cluster->*accessor)(1,1) +
      2.*cosphi*sinphi*(cluster->*accessor)(0,1);
  }
}

TrkrClusterv3::TrkrClusterv3()
  : m_cluskey(TrkrDefs::CLUSKEYMAX)
  , m_subsurfkey(TrkrDefs::SUBSURFKEYMAX)
  , m_adc(0xFFFFFFFF)
{
  for (int i = 0; i < 6; ++i)
    { m_size[i] = NAN; }

  for (int i = 0; i < 2; i++)
    {
      m_local[i] = NAN;
      for(int j = 0; j < 2; j ++)
	{
	  m_actsLocalErr[i][j] = NAN;
	}
    }
}

void TrkrClusterv3::identify(std::ostream& os) const
{
  os << "---TrkrClusterv3--------------------" << std::endl;
  os << "clusid: " << getClusKey() << std::dec << std::endl;

  os << " (rphi,z) =  (" << getLocalX();
  os << ", " << getLocalY() << ") cm ";

  os << " adc = " << getAdc() << std::endl;

  os << "Size " << std::endl;
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
  os << getError(0, 1) << " ) " << std::endl;
  os << "         ( ";
  os << getError(1, 0) << ", ";
  os << getError(1, 1) << " ) " << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int TrkrClusterv3::isValid() const
{
  if (m_cluskey == TrkrDefs::CLUSKEYMAX) { return 0; }
  for (int i = 0; i < 2; ++i)
  {
    if (std::isnan(getPosition(i))) { return 0; }
    for (int j = 0; j < 2; ++j)
      { if (std::isnan(getActsLocalError(i,j))) { return 0; }}
  }
  if (m_adc == 0xFFFFFFFF) { return 0; }
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (std::isnan(getSize(i, j))) { return 0; }
    }
  }

  return 1;
}

void TrkrClusterv3::setSize(unsigned int i, unsigned int j, float value)
{
  m_size[covarIndex(i, j)] = value;
  return;
}

float TrkrClusterv3::getSize(unsigned int i, unsigned int j) const
{ return m_size[covarIndex(i, j)]; }


float TrkrClusterv3::getPhiSize() const
{ return 2*std::sqrt(rotate<&TrkrClusterv3::getSize>(this)); }

float TrkrClusterv3::getZSize() const
{ return 2.*sqrt(getSize(2, 2)); }

float TrkrClusterv3::getRPhiError() const
{ return std::sqrt(m_actsLocalErr[0][0]); }

float TrkrClusterv3::getZError() const
{ return std::sqrt(m_actsLocalErr[1][1]); }

void TrkrClusterv3::setActsLocalError(unsigned int i, unsigned int j,
				      float value)
{
  m_actsLocalErr[i][j] = value;
}
