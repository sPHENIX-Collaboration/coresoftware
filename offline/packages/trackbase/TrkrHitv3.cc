#include "TrkrHitv3.h"
#include <climits>

TrkrHitv3::TrkrHitv3()
{
}

  // these set and get the energy before digitization
void TrkrHitv3::addEnergy(const double edep) 
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

double TrkrHitv3::getEnergy() 
{
  return ((double) m_adc)  / TrkrDefs::EdepScaleFactor;
}

void TrkrHitv3::setAdc(const unsigned int adc)
 {
   if(adc > USHRT_MAX)
     m_adc = USHRT_MAX;
   else
     m_adc = (unsigned short) adc;
 }

unsigned int TrkrHitv3::getAdc() { 
    return (unsigned int) m_adc;
  }

void TrkrHitv3::setCrossing(const short int crossing) { 
    m_crossing = crossing;
  }

short int TrkrHitv3::getCrossing() { 
    return m_crossing;
  }
