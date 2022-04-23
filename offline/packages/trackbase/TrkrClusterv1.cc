/**
 * @file trackbase/TrkrClusterv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrClusterv1
 */
#include "TrkrClusterv1.h"

#include <cmath>
#include <utility>          // for swap

namespace
{

  // square convenience function
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  // get unique index in cov. matrix array from i and j
  inline unsigned int covarIndex(unsigned int i, unsigned int j)
  {
    if (i > j) std::swap(i, j);
    return i + 1 + (j + 1) * (j) / 2 - 1;
  }

 // rotate size or covariant matrix to polar coordinates and return the phi component
 template<float (TrkrClusterv1::*accessor)(unsigned int, unsigned int) const>
    float rotate( const TrkrClusterv1* cluster )
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

TrkrClusterv1::TrkrClusterv1()
  : m_cluskey(TrkrDefs::CLUSKEYMAX)
  , m_isGlobal(true)
  , m_adc(0xFFFFFFFF)
{
  for (int i = 0; i < 3; ++i) m_pos[i] = NAN;

  for (int j = 0; j < 3; ++j)
  {
    for (int i = 0; i < 3; ++i)
    {
      setSize(i, j, NAN);
      setError(i, j, NAN);
    }
  }
}

void TrkrClusterv1::identify(std::ostream& os) const
{
  os << "---TrkrClusterv1--------------------" << std::endl;

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
    if (std::isnan(getPosition(i))) return 0;
  }
  if (m_adc == 0xFFFFFFFF) return 0;
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (std::isnan(getSize(i, j))) return 0;
      if (std::isnan(getError(i, j))) return 0;
    }
  }

  return 1;
}

void TrkrClusterv1::CopyFrom( const TrkrCluster& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
 
  // parent class method
  TrkrCluster::CopyFrom( source );

  setX( source.getX() );
  setY( source.getY() );
  setZ( source.getZ() );
  m_isGlobal = source.isGlobal();
  setAdc( source.getAdc() );

  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
  {
    setSize(i, j, source.getSize(i, j));
    setError(i, j, source.getError(i, j));
  }
}

void TrkrClusterv1::setSize(unsigned int i, unsigned int j, float value)
{
  m_size[covarIndex(i, j)] = value;
  return;
}

float TrkrClusterv1::getSize(unsigned int i, unsigned int j) const
{ return m_size[covarIndex(i, j)]; }

void TrkrClusterv1::setError(unsigned int i, unsigned int j, float value)
{
  m_err[covarIndex(i, j)] = value;
  return;
}

float TrkrClusterv1::getError(unsigned int i, unsigned int j) const
{ return m_err[covarIndex(i, j)]; }

float TrkrClusterv1::getPhiSize() const
{ return 2*std::sqrt(rotate<&TrkrClusterv1::getSize>(this)); }

float TrkrClusterv1::getZSize() const
{ return 2.*sqrt(getSize(2, 2)); }

float TrkrClusterv1::getPhiError() const
{
  const float rad = std::sqrt(square(m_pos[0])+square(m_pos[1]));
  if (rad > 0) return getRPhiError() / rad;
  return 0;
}

float TrkrClusterv1::getRPhiError() const
{ return std::sqrt(rotate<&TrkrClusterv1::getError>( this )); }

float TrkrClusterv1::getZError() const
{ return std::sqrt(getError(2, 2)); }
