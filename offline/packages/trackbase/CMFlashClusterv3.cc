/**
 * @file trackbase/CMFlashClusterv3.cc
 * @author Ben Kimelman
 * @date March 2023
 * @brief Implementation of CMFlashClusterv3
 */
#include "CMFlashClusterv3.h"

#include <cmath>
#include <utility>  // for swap

void CMFlashClusterv3::identify(std::ostream& os) const
{
  os << "---CMFlashClusterv3--------------------" << std::endl;

  os << " (x,y,z) =  (" << m_pos[0];
  os << ", " << m_pos[1] << ", ";
  os << m_pos[2] << ") cm";

  os << " adc = " << getAdc() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int CMFlashClusterv3::isValid() const
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

  if (std::isnan(getX1()))
  {
    return 0;
  }
  if (std::isnan(getY1()))
  {
    return 0;
  }
  if (std::isnan(getZ1()))
  {
    return 0;
  }

  if (std::isnan(getX2()))
  {
    return 0;
  }
  if (std::isnan(getY2()))
  {
    return 0;
  }
  if (std::isnan(getZ2()))
  {
    return 0;
  }

  if (m_adc == 0xFFFFFFFF)
  {
    return 0;
  }
  if (m_adc1 == 0xFFFFFFFF)
  {
    return 0;
  }
  if (m_adc2 == 0xFFFFFFFF)
  {
    return 0;
  }

  return 1;
}

void CMFlashClusterv3::CopyFrom(const CMFlashCluster& source)
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

  setX1(source.getX1());
  setY1(source.getY1());
  setZ1(source.getZ1());

  setX2(source.getX2());
  setY2(source.getY2());
  setZ2(source.getZ2());

  setLayer1(source.getLayer1());
  setLayer2(source.getLayer2());

  setAdc(source.getAdc());
  setAdc1(source.getAdc1());
  setAdc2(source.getAdc2());
  setIsRGap(source.getIsRGap());
  setIsPhiGap(source.getIsPhiGap());
}
