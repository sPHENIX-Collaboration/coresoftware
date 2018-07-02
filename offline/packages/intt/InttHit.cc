/**
 * @file intt/InttHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Intt hit object
 */
#include "InttHit.h"
#include "InttDefs.h"

InttHit::InttHit()
  : TrkrHit()
  , m_adc(0)
{
}

void 
InttHit::identify(std::ostream& os) const
{
  os << "InttHit with adc:" << getAdc()
     << std::endl;
}

void 
InttHit::Reset()
{
  TrkrHit::Reset();
}

int 
InttHit::isValid() const
{
  // valid if the key is not equal to the default value
  return m_adc != 0;
}

