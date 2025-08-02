/**
 * @file trackbase/LaserClusterv2.cc
 * @author Ben Kimelman
 * @date July 2025
 * @brief Implementation of LaserClusterv2
 */
#include "LaserClusterv2.h"

#include <cmath>
#include <utility>          // for swap

void LaserClusterv2::identify(std::ostream& os) const
{
  os << "---LaserClusterv2--------------------" << std::endl;

  os << " fit? " << m_fitMode;
  os << " (x,y,z) =  (" << m_pos[0];
  os << ", " << m_pos[1] << ", ";
  os << m_pos[2] << ") cm";

  os << " adc = " << getAdc() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int LaserClusterv2::isValid() const
{
  if(std::isnan(getX()))
    {
      return 0;
    }
  if(std::isnan(getY()))
    {
      return 0;
    }
  if(std::isnan(getZ()))
    {
      return 0;
    }

  if (m_adc == 0xFFFFFFFF)
    {
      return 0;
    }

  return 1;
}

void LaserClusterv2::CopyFrom( const LaserCluster& source )
{
  // do nothing if copying onto oneself
  if( this == &source )
    {
      return;
    }
 
  // parent class method
  LaserCluster::CopyFrom( source );

  setX( source.getX() );
  setY( source.getY() );
  setZ( source.getZ() );
  setFitMode( source.getFitMode() );
  setAdc( source.getAdc() );
  setLayer( source.getLayer() );
  setIPhi( source.getIPhi() );
  setIT( source.getIT() );
  setNhits( source.getNhits() );
  setNLayers( source.getNLayers() );
  setNIPhi( source.getNIPhi() );
  setNIT( source.getNIT() );
  setSDLayer( source.getSDLayer() );
  setSDIPhi( source.getSDIPhi() );
  setSDIT( source.getSDIT() );
  setSDWeightedLayer( source.getSDWeightedLayer() );
  setSDWeightedIPhi( source.getSDWeightedIPhi() );
  setSDWeightedIT( source.getSDWeightedIT() );


  for(int i=0; i<(int)source.getNhits(); i++){
    addHit();
    setHitLayer(i, source.getHitLayer(i));
    setHitIPhi(i, source.getHitIPhi(i));
    setHitIT(i, source.getHitIT(i));

    setHitX(i, source.getHitX(i));
    setHitY(i, source.getHitY(i));
    setHitZ(i, source.getHitZ(i));
    setHitAdc(i, source.getHitAdc(i));
  }


}

