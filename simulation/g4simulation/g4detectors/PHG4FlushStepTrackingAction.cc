#include "PHG4FlushStepTrackingAction.h"

#include <g4main/PHG4SteppingAction.h>

#include <Geant4/G4Track.hh>

#include <iostream>

using namespace std;

PHG4FlushStepTrackingAction::PHG4FlushStepTrackingAction(PHG4SteppingAction *stepact):
steppingAction_(stepact)
{}

void
PHG4FlushStepTrackingAction::PostUserTrackingAction(const G4Track* trk)
{
  if (steppingAction_)
    {
      steppingAction_->flush_cached_values();
    }
  return;
}
