#include "PHG4ConeRegionSteppingAction.h"

#include "PHG4ConeDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>              // for G4StepPoint
#include <Geant4/G4StepStatus.hh>             // for fGeomBoundary, fUndefined
#include <Geant4/G4SystemOfUnits.hh>     // for cm
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>        // for G4TouchableHandle
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <cstddef>                           // for NULL
#include <iostream>
#include <string>                             // for operator+, operator<<

class G4VPhysicalVolume;

using namespace std;
//____________________________________________________________________________..
PHG4ConeRegionSteppingAction::PHG4ConeRegionSteppingAction( PHG4ConeDetector* detector ):
  detector_( detector ), hits_(nullptr), hit(nullptr)
{}

//____________________________________________________________________________..
void PHG4ConeRegionSteppingAction::UserSteppingAction( const G4Step* aStep)
{

  // get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;

  const G4Track* aTrack = aStep->GetTrack();

  // make sure we are in a volume
  if ( detector_->IsInConeActive(volume) )
    {
      int layer_id = 0;
      G4StepPoint * prePoint = aStep->GetPreStepPoint();
      G4StepPoint * postPoint = aStep->GetPostStepPoint();
       cout << "track id " << aTrack->GetTrackID() << endl;
       cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
       cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
      switch (prePoint->GetStepStatus())
        {
        case fGeomBoundary:
        case fUndefined:
          hit = new PHG4Hitv1();
          //here we set the entrance values in cm
          hit->set_x( 0, prePoint->GetPosition().x() / cm);
          hit->set_y( 0, prePoint->GetPosition().y() / cm );
          hit->set_z( 0, prePoint->GetPosition().z() / cm );
	  // time in ns
          hit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );
          //set the track ID
	  {
	    hit->set_trkid(aTrack->GetTrackID());
	    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	      {
		if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		  {
		    hit->set_trkid(pp->GetUserTrackId());
		    hit->set_shower_id(pp->GetShower()->get_id());
		  }
	      }
	  }
	  
          //set the initial energy deposit
          hit->set_edep(0);

          // Now add the hit
          hits_->AddHit(layer_id, hit);

	  {
	    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	      {
		if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		  {
		    pp->GetShower()->add_g4hit_id(hits_->GetID(),hit->get_hit_id());
		  }
	      }
	  }

          break;
        default:
          break;
        }
      // here we just update the exit values, it will be overwritten
      // for every step until we leave the volume or the particle
      // ceases to exist
      hit->set_x( 1, postPoint->GetPosition().x() / cm );
      hit->set_y( 1, postPoint->GetPosition().y() / cm );
      hit->set_z( 1, postPoint->GetPosition().z() / cm );

      hit->set_t( 1, postPoint->GetGlobalTime() / nanosecond );
      //sum up the energy to get total deposited
      hit->set_edep(hit->get_edep() + edep);

      //    hit->identify();
      // return true to indicate the hit was used
      return;

    }
  else
    {
      return;
    }
}

//____________________________________________________________________________..
void PHG4ConeRegionSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
{

  string hitnodename = "G4HIT_" + detector_->GetName();
  //now look for the map and grab a pointer to it.
  hits_ =  findNode::getClass<PHG4HitContainer>( topNode , hitnodename.c_str() );

  // if we do not find the node we need to make it.
  if ( ! hits_ )
    { std::cout << "PHG4ConeRegionSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
 }

}
