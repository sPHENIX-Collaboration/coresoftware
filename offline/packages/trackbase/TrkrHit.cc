#include "TrkrHit.h"
#include <climits>

TrkrHit::TrkrHit()
  :  m_adc(0)
{
}

  // these set and get the energy before digitization
void TrkrHit::addEnergy(const double edep) 
{
  double max_adc = (double) USHRT_MAX;

  // overflow occurs if (sum of the existing + new ADC values)   > USHRT_MAX
  // where new ADC values are (edep * TrkrDefs::EdepScaleFactor)
  double ein = edep* TrkrDefs::EdepScaleFactor;
  if( (double) m_adc + ein > max_adc)  
    {
    m_adc = USHRT_MAX;
    //std::cout << " warning: exceeded range of unsigned short: ein " << ein << " m_adc " << m_adc << " max_adc " << max_adc << std::endl;
    }
  else
    {
      m_adc += (unsigned short) (ein); 
      //std::cout << " edep " << edep << " ein " << ein << " m_adc " << m_adc << " max_adc " << max_adc << std::endl;  
    }
}

double TrkrHit::getEnergy() 
{
  //std::cout << " return m_adc = " << m_adc << " energy " <<  ((double) m_adc)  / TrkrDefs::EdepScaleFactor << std::endl;
  return ((double) m_adc)  / TrkrDefs::EdepScaleFactor;
}
