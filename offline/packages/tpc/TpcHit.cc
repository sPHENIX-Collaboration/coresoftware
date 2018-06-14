/**
 * @file tpc/TpcHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Tpc hit object
 */
#include "TpcHit.h"
#include "TpcDefs.h"

TpcHit::TpcHit()
  : TrkrHit()
{
}

void 
TpcHit::identify(std::ostream& os) const
{
  os << "TpcHit with key:" << getKey() 
     << " and adc:" << getAdc() 
     << std::endl;
}

void 
TpcHit::Reset()
{
  TrkrHit::Reset();
}

int 
TpcHit::isValid() const
{
  // valid if the key is not equal to the default value
  return getKey() != TrkrDefs::HITKEYMAX;
}

