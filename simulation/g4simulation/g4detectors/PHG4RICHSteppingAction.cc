// $$Id: PHG4RICHSteppingAction.cc,v 1.5 2015/01/07 21:54:33 pinkenbu Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov> and Nils Feege <nils.feege@stonybrook.edu>
 * \version $$Revision: 1.5 $$
 * \date $$Date: 2015/01/07 21:54:33 $$
 */

#include "PHG4RICHSteppingAction.h"
#include "PHG4RICHDetector.h"
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4SDManager.hh"

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>

using namespace std;

PHG4RICHSteppingAction::PHG4RICHSteppingAction(PHG4RICHDetector* detector) :
  detector_(detector),
  hits_(NULL),
  hit(NULL)
{

}

void
PHG4RICHSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* theTrack = aStep->GetTrack();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4OpBoundaryProcess* boundary=NULL;

  /* find the boundary process only once */
  if(!boundary){
    G4ProcessManager* pm  = aStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

  if(!thePostPV){ //out of world
    return;
  }

  /* Optical photon only */
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){

    /* Was the photon absorbed by the absorption process? */
    if(thePostPoint->GetProcessDefinedStep()->GetProcessName() == "OpAbsorption")
      {
      }
    assert(boundary != NULL);
    boundaryStatus=boundary->GetStatus();

    /*Check to see if the partcile was actually at a boundary
      Otherwise the boundary status may not be valid
      Prior to Geant4.6.0-p1 this would not have been enough to check  */
    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      if(fExpectedNextStatus==StepTooSmall){
        if(boundaryStatus!=StepTooSmall){
          G4ExceptionDescription ed;
          ed << "EicRichGemTbSteppingAction::UserSteppingAction(): "
             << "No reallocation step after reflection!"
             << G4endl;
          G4Exception("EicRichGemTbSteppingAction::UserSteppingAction()", "EicRichGemTbExpl01",
                      FatalException,ed,
                      "Something is wrong with the surface normal or geometry");
        }
      }
      fExpectedNextStatus=Undefined;
      switch(boundaryStatus){
      case Absorption:
        break;
      case Detection: /*Note, this assumes that the volume causing detection
                        is the photocathode because it is the only one with
                        non-zero efficiency */
        {
          /* make sure the photon actually did hit the GEM stack */
	  if ( thePostPV->GetName() == "RICHHBDGEMFrontCu1Physical" )
            {
              MakeHit( aStep );
            }
          break;
        }
      case FresnelReflection:
      case TotalInternalReflection:
      case LambertianReflection:
      case LobeReflection:
      case SpikeReflection:
      case BackScattering:
        fExpectedNextStatus=StepTooSmall;
        break;
      default:
        break;
      }
    }
  }

  return;
}

bool PHG4RICHSteppingAction::MakeHit(const G4Step* aStep){

  const G4Track* aTrack = aStep->GetTrack();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();

  int layer_id = 0;

  hit = new PHG4Hitv1();

  //here we set the entrance values in cm
  hit->set_x( 0, postPoint->GetPosition().x() / cm);
  hit->set_y( 0, postPoint->GetPosition().y() / cm );
  hit->set_z( 0, postPoint->GetPosition().z() / cm );

  // time in ns
  hit->set_t( 0, postPoint->GetGlobalTime() / nanosecond );

  //same for exit values (photons absorbed/detected at boundary to post step volume)
  hit->set_x( 1, postPoint->GetPosition().x() / cm );
  hit->set_y( 1, postPoint->GetPosition().y() / cm );
  hit->set_z( 1, postPoint->GetPosition().z() / cm );

  hit->set_t( 1, postPoint->GetGlobalTime() / nanosecond );

  //set the track ID
  {
    hit->set_trkid(aTrack->GetTrackID());
    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
      {
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	  {
	    hit->set_trkid(pp->GetUserTrackId());
	  }
      }
  }

  // set optical photon energy deposition to 0
  hit->set_edep(0);

  // Now add the hit
  hits_->AddHit(layer_id, hit);

  // return true to indicate the hit was used
  return true;
}

void
PHG4RICHSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{

  //now look for the map and grab a pointer to it.
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_RICH");

  // if we do not find the node we need to make it.
  if (!hits_)
    {
      std::cout
        << "PHG4RICHSteppingAction::SetTopNode - unable to find G4HIT_RICH"
        << std::endl;
    }

}

