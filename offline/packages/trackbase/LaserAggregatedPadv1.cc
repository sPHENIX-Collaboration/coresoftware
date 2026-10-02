/**
 * @file trackbase/LaserAggregatedPadv1.cc
 * @author Ben Kimelman
 * @date July 2026
 * @brief Implementation of LaserAggregatedPadv1
 */
#include "LaserAggregatedPadv1.h"

#include <cmath>
#include <utility>          // for swap

void LaserAggregatedPadv1::identify(std::ostream& os) const
{
  os << "---LaserAggregatedPadv1--------------------" << std::endl;
  os << " adc = " << getAdc() << "   nHits: " << getNHits() << std::endl;

  os << std::endl;
  os << "-----------------------------------------------" << std::endl;

  return;
}

void LaserAggregatedPadv1::Reset()
{
  m_adc = 0.0;
  m_nHits = 0.0;
}

int LaserAggregatedPadv1::isValid() const
{
  if(getNHits() == 0)
  {
    return 0;
  }

  return 1;
}

void LaserAggregatedPadv1::CopyFrom( const LaserAggregatedPad& source )
{
  if( this == &source ) return;
  LaserAggregatedPad::CopyFrom( source );
  m_adc   = source.getAdc();
  m_nHits = source.getNHits();
}

