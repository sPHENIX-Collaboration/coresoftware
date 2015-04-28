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
//   if ( VERBOSE>0 )
//     {
//       G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
//       std::cout << "PHG4TruthSteppingAction::UserSteppingAction: Stepping in volume " << volume->GetName()
// 		<< " = " << volume;
//       if ( volume->IsParameterised() ) std::cout << " copyNo " << volume->GetCopyNo();
//       std::cout << " " << (hitWasUsed?"hit was used":"hit was not used") << std::endl;
//     }

//   bool addid = false;

//   const G4Track* track=step->GetTrack();
//   if ( ! track ) return false;

//   //  if( !hitWasUsed ) return false;

//   // hitWasUsed being true means that the hit corresponds to at least
//   // one of the previously defined subsystem, in which case one adds its
//   // corresponding track id to the eventAction
//   // We have to check both hitused and the user data stored with the track, since 
//   // region stepping routines will be called outside the master stepping action.

//   if ( hitWasUsed ) addid = true;
//   else if ( PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()) )
//     {
//       if ( p->GetWanted() ) addid = true;
//     }

//   if ( ! addid ) return false;

//   if( track )
//   {
//     eventAction_->AddTrackidToWritelist( track->GetTrackID() );
//     return true;

//   } else return false;
  return false;

}
