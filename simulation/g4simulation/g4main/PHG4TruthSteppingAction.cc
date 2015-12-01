#include "PHG4TruthSteppingAction.h"
#include "PHG4TrackUserInfoV1.h"
#include "PHG4TruthEventAction.h"

#include <Geant4/G4Step.hh>
#include <Geant4/G4Track.hh>

const int VERBOSE = 0;

//________________________________________________________
PHG4TruthSteppingAction::PHG4TruthSteppingAction( PHG4TruthEventAction* eventAction ):
  eventAction_( eventAction )
  {}

//________________________________________________________
bool PHG4TruthSteppingAction::UserSteppingAction( const G4Step* step, bool hitWasUsed )
{
  return false;
}
