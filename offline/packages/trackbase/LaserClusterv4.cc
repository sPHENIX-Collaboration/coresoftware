/**
 * @file trackbase/LaserClusterv4.cc
 * @author Ben Kimelman
 * @date July 2026
 * @brief Implementation of LaserClusterv4
 */
#include "LaserClusterv4.h"

#include <cmath>
#include <utility>          // for swap

void LaserClusterv4::identify(std::ostream& os) const
{
  os << "---LaserClusterv4--------------------" << std::endl;

  os << " " << m_hits.size() << " hits";
  os << " fit? " << m_fitMode;
  os << " adc = " << getAdc() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int LaserClusterv4::isValid() const
{
  if(getNhits() == 0)
  {
    return 0;
  }

  return 1;
}

unsigned int LaserClusterv4::getAdc() const
{
  unsigned int adc = 0;
  for(const auto &LCHI : m_hits)
  {
    adc += (unsigned int) LCHI.adc;
  }
  return adc;
}

void LaserClusterv4::CopyFrom( const LaserCluster& source )
{
  // do nothing if copying onto oneself
  if( this == &source )
    {
      return;
    }
 
  // parent class method
  LaserCluster::CopyFrom( source );
  m_hits.clear();
  setFitMode(source.getFitMode());
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

