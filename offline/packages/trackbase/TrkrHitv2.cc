#include "TrkrHitv2.h"
#include <climits>

TrkrHitv2::TrkrHitv2()
{
}

  // these set and get the energy before digitization
void TrkrHitv2::addEnergy(const double edep) 
{
  double max_adc = (double) USHRT_MAX;

  // overflow occurs if (sum of the existing + new ADC values)   > USHRT_MAX
  double ein = edep* TrkrDefs::EdepScaleFactor;
  if( (double) m_adc + ein > max_adc)  
    {
    m_adc = USHRT_MAX;
    }
  else
    {
      m_adc += (unsigned short) (ein); 
    }
}

double TrkrHitv2::getEnergy() 
{
  return ((double) m_adc)  / TrkrDefs::EdepScaleFactor;
}

void TrkrHitv2::setAdc(const unsigned int adc)
 {
   if(adc > USHRT_MAX)
     m_adc = USHRT_MAX;
   else
     m_adc = (unsigned short) adc;
 }

unsigned int TrkrHitv2::getAdc() { 
    return (unsigned int) m_adc;
  }
