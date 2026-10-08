/**
 * @file trackbase/LaserClusterv5.cc
 * @author Ben Kimelman
 * @date Oct 2026
 * @brief Implementation of LaserClusterv5
 */
#include "LaserClusterv5.h"

#include <cmath>
#include <utility>          // for swap

void LaserClusterv5::identify(std::ostream& os) const
{
  os << "---LaserClusterv5--------------------" << std::endl;

  os << " " << m_hits.size() << " hits";
  os << " " << m_hitsDouble.size() << " double stored adc hits";
  os << " fit? " << m_fitMode;
  os << " adc = " << getAdc();
  os << " double stored adc = " << getAdcDouble();
  os << " isLamination? " << getIsLamination();
  os << " truthIndex = " << getTruthIndex() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

int LaserClusterv5::isValid() const
{
  if(getNhits() == 0 && getNhitsDouble() == 0)
  {
    return 0;
  }

  return 1;
}

unsigned int LaserClusterv5::getAdc() const
{
  unsigned int adc = 0;
  for(const auto &LCHI : m_hits)
  {
    adc += (unsigned int) LCHI.adc;
  }
  return adc;
}

double LaserClusterv5::getAdcDouble() const
{
  double adc = 0;
  for(const auto &LCHI : m_hitsDouble)
  {
    adc += (double) LCHI.adc;
  }
  return adc;
}

void LaserClusterv5::CopyFrom( const LaserCluster& source )
{
  // do nothing if copying onto oneself
  if( this == &source )
    {
      return;
    }
 
  // parent class method
  LaserCluster::CopyFrom( source );
  m_hits.clear();
  m_hitsDouble.clear();
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

  setTruthIndex( source.getTruthIndex() );
  setIsLamination( source.getIsLamination() );


  for(int i=0; i<(int)source.getNhits(); i++){
    LaserClusterHitInfo LCHI = source.getHit(i);
    addHit(LCHI.hitsetkey, LCHI.hitkey, LCHI.adc);
  }

  for(int i=0; i<(int)source.getNhitsDouble(); i++){
    LaserClusterHitInfoDouble LCHI = source.getHitDouble(i);
    addHitDouble(LCHI.hitsetkey, LCHI.hitkey, LCHI.adc);
  }
}

