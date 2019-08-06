/**
 * @file outertracker/OuterTrackerHit.cc
 * @author A Frawley
 * @date June 2018
 * @brief Implementation of OuterTracker hit object
 */
#include "OuterTrackerHit.h"

#include <trackbase/TrkrHit.h>  // for TrkrHit

#include <iosfwd>

OuterTrackerHit::OuterTrackerHit()
  : TrkrHit()
{
}

void OuterTrackerHit::identify(std::ostream& os) const
{
  //os << "OuterTrackerHit with energy " << getEnergy() << " adc:" << getAdc()
  // << std::endl;
}

void OuterTrackerHit::Reset()
{
  TrkrHit::Reset();
}

int OuterTrackerHit::isValid() const
{
  // valid if the key is not equal to the default value
  return 0;
}
