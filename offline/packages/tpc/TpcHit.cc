/**
 * @file tpc/TpcHit.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of Tpc hit object
 */
#include "TpcHit.h"

#include <trackbase/TrkrHit.h>  // for TrkrHit

TpcHit::TpcHit()
  : TrkrHit()
{
}

void TpcHit::identify(std::ostream& os) const
{

}

void TpcHit::Reset()
{
  TrkrHit::Reset();
}

int TpcHit::isValid() const
{
  // valid if the adc is not equal to the default value
  return 0;
}
