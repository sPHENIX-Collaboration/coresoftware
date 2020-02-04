/**
 * @file intt/InttHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Intt hit object
 */
#include "InttHit.h"

#include <trackbase/TrkrHit.h>  // for TrkrHit

InttHit::InttHit()
  : TrkrHit()
{
}

void InttHit::identify(std::ostream& os) const
{
  //os << "InttHit with energy " << getEnergy() << " adc:" << getAdc()
  // << std::endl;
}

void InttHit::Reset()
{
  TrkrHit::Reset();
}

int InttHit::isValid() const
{
  // valid if the key is not equal to the default value
  return 0;
}
