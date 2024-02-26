/**
 * @file trackbase/TrkrClusterv5.cc
 * @author J. Osborn
 * @date October 2021
 * @brief Implementation of TrkrClusterv5
 */
#include "TrkrClusterv5.h"

#include <cmath>
#include <utility>  // for swap

namespace
{
  // square convenience function
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

TrkrClusterv5::TrkrClusterv5()
  : m_subsurfkey(TrkrDefs::SUBSURFKEYMAX)
  , m_phierr(0)
  , m_zerr(0)
  , m_adc(0)
  , m_maxadc(0)
  , m_phisize(0)
  , m_zsize(0)
  , m_overlap(0)
  , m_edge(0)
{
  for (float& i : m_local)
  {
    i = NAN;
  }
}

void TrkrClusterv5::identify(std::ostream& os) const
{
  os << "---TrkrClusterv5--------------------" << std::endl;

  os << " (rphi,z) =  (" << getLocalX();
  os << ", " << getLocalY() << ") cm ";

  os << " valid = " << isValid() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int TrkrClusterv5::isValid() const
{
  for (int i = 0; i < 2; ++i)
  {
    if (std::isnan(getPosition(i)))
    {
      return 0;
    }
  }
  if (m_adc == 0xFFFF)
  {
    return 0;
  }

  return 1;
}

void TrkrClusterv5::CopyFrom(const TrkrCluster& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  TrkrCluster::CopyFrom(source);

  setLocalX(source.getLocalX());
  setLocalY(source.getLocalY());
  setSubSurfKey(source.getSubSurfKey());
  setAdc(source.getAdc());
  setMaxAdc(source.getMaxAdc());
  setPhiError(source.getPhiError());
  setZError(source.getZError());
  setPhiSize(source.getPhiSize());
  setZSize(source.getZSize());
  setOverlap(source.getOverlap());
  setEdge(source.getEdge());
}
