/**
 * @file trackbase/LaserClusterv3.cc
 * @author Ben Kimelman
 * @date July 2026
 * @brief Implementation of LaserClusterv3
 */
#include "LaserClusterv3.h"

#include <cmath>
#include <utility>          // for swap

void LaserClusterv3::identify(std::ostream& os) const
{
  os << "---LaserClusterv3--------------------" << std::endl;

  os << " " << m_hits.size() << " hits";
  os << " fit? " << m_fitMode;
  os << " (layer, iphi, it) =  (" << m_posHardware[0];
  os << ", " << m_posHardware[1] << ", ";
  os << m_posHardware[2] << ")";
  os << " adc = " << getAdc() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int LaserClusterv3::isValid() const
{
  if(getNhits() == 0)
  {
    return 0;
  }

  return 1;
}

unsigned int LaserClusterv3::getAdc() const
{
  unsigned int adc = 0;
  for(const auto &LCHI : m_hits)
  {
    adc += (unsigned int) LCHI.adc;
  }
  return adc;
}

void LaserClusterv3::CopyFrom( const LaserCluster& source )
{
  // do nothing if copying onto oneself
  if( this == &source )
    {
      return;
    }
 
  // parent class method
  LaserCluster::CopyFrom( source );
  setLayerInt( source.getLayerInt() );
  setIPhiInt( source.getIPhiInt() );
  setITInt( source.getITInt() );
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
    LaserClusterHitInfo LCHI = source.getHit(i);
    addHit(LCHI.hitsetkey, LCHI.hitkey, LCHI.adc);
  }
}

