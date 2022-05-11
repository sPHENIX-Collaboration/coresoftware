/**
 * @file trackbase/CMFlashClusterv1.cc
 * @author Tony Frawley
 * @date January 2022
 * @brief Implementation of CMFlashClusterv1
 */
#include "CMFlashClusterv1.h"

#include <cmath>
#include <utility>          // for swap

void CMFlashClusterv1::identify(std::ostream& os) const
{
  os << "---CMFlashClusterv1--------------------" << std::endl;

  os << " (x,y,z) =  (" << m_pos[0];
  os << ", " << m_pos[1] << ", ";
  os << m_pos[2] << ") cm";

  os << " adc = " << getAdc() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int CMFlashClusterv1::isValid() const
{
  if(std::isnan(getX())) return 0;
  if(std::isnan(getY())) return 0;
  if(std::isnan(getZ())) return 0;

  if (m_adc == 0xFFFFFFFF) return 0;

  return 1;
}

void CMFlashClusterv1::CopyFrom( const CMFlashCluster& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
 
  // parent class method
  CMFlashCluster::CopyFrom( source );

  setX( source.getX() );
  setY( source.getY() );
  setZ( source.getZ() );
  setAdc( source.getAdc() );

}

