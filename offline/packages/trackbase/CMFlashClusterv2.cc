/**
 * @file trackbase/CMFlashClusterv2.cc
 * @author Ben Kimelman
 * @date March 2023
 * @brief Implementation of CMFlashClusterv2
 */
#include "CMFlashClusterv2.h"

#include <cmath>
#include <utility>  // for swap

void CMFlashClusterv2::identify(std::ostream& os) const
{
  os << "---CMFlashClusterv2--------------------" << std::endl;

  os << " (x,y,z) =  (" << m_pos[0];
  os << ", " << m_pos[1] << ", ";
  os << m_pos[2] << ") cm";

  os << " adc = " << getAdc() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int CMFlashClusterv2::isValid() const
{
  if (std::isnan(getX()))
  {
    return 0;
  }
  if (std::isnan(getY()))
  {
    return 0;
  }
  if (std::isnan(getZ()))
  {
    return 0;
  }

  if (m_adc == 0xFFFFFFFF)
  {
    return 0;
  }

  return 1;
}

void CMFlashClusterv2::CopyFrom(const CMFlashCluster& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  CMFlashCluster::CopyFrom(source);

  setX(source.getX());
  setY(source.getY());
  setZ(source.getZ());
  setAdc(source.getAdc());
  setIsRGap(source.getIsRGap());
  setIsPhiGap(source.getIsPhiGap());
}
