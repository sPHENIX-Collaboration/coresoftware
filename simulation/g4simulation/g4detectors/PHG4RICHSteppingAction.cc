/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov> and Nils Feege <nils.feege@stonybrook.edu>
 * \version $$Revision: 1.5 $$
 * \date $$Date: 2015/01/07 21:54:33 $$
 */

#include "PHG4RICHSteppingAction.h"

#include "PHG4RICHDetector.h"
#include "ePHENIXRICHConstruction.h"          // for ePHENIXRICHConstruction

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4ExceptionSeverity.hh>      // for FatalException
#include <Geant4/G4OpticalPhoton.hh>          // for G4OpticalPhoton
#include <Geant4/G4ParticleDefinition.hh>     // for G4ParticleDefinition
#include <Geant4/G4ProcessManager.hh>
#include <Geant4/G4ProcessVector.hh>          // for G4ProcessVector
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>              // for G4StepPoint
#include <Geant4/G4StepStatus.hh>             // for fGeomBoundary
#include <Geant4/G4String.hh>                 // for G4String
#include <Geant4/G4SystemOfUnits.hh>          // for cm, nanosecond, GeV
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4Types.hh>                  // for G4int, G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VProcess.hh>               // for G4VProcess
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation
#include <Geant4/G4ios.hh>                    // for G4endl
#include <Geant4/globals.hh>                  // for G4Exception, G4Exceptio...

#include <cassert>                           // for assert
#include <iostream>                           // for operator<<, basic_ostream

using namespace std;

PHG4RICHSteppingAction::PHG4RICHSteppingAction(PHG4RICHDetector* detector)
  : detector_(detector)
  , hits_(nullptr)
  , hit(nullptr)
  , fExpectedNextStatus(Undefined)
{
}

void PHG4RICHSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* theTrack = aStep->GetTrack();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus = Undefined;
  static G4OpBoundaryProcess* boundary = nullptr;

  /* find the boundary process only once */
  if (!boundary)
  {
    G4ProcessManager* pm = aStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for (i = 0; i < nprocesses; i++)
    {
      if ((*pv)[i]->GetProcessName() == "OpBoundary")
      {
        boundary = (G4OpBoundaryProcess*) (*pv)[i];
        break;
      }
    }
  }

  if (!thePostPV)
  {  //out of world
    return;
  }

  /* Optical photon only */
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
  {
    /* Was the photon absorbed by the absorption process? */
    if (thePostPoint->GetProcessDefinedStep()->GetProcessName() == "OpAbsorption")
    {
    }
    assert(boundary != nullptr);
    boundaryStatus = boundary->GetStatus();

    /*Check to see if the partcile was actually at a boundary
      Otherwise the boundary status may not be valid
      Prior to Geant4.6.0-p1 this would not have been enough to check  */
    if (thePostPoint->GetStepStatus() == fGeomBoundary)
    {
      if (fExpectedNextStatus == StepTooSmall)
      {
        if (boundaryStatus != StepTooSmall)
        {
          G4ExceptionDescription ed;
          ed << "EicRichGemTbSteppingAction::UserSteppingAction(): "
             << "No reallocation step after reflection!"
             << G4endl;
          G4Exception("EicRichGemTbSteppingAction::UserSteppingAction()", "EicRichGemTbExpl01",
                      FatalException, ed,
                      "Something is wrong with the surface normal or geometry");
        }
      }
      fExpectedNextStatus = Undefined;
      switch (boundaryStatus)
      {
      case Absorption:
        break;
      case Detection: /*Note, this assumes that the volume causing detection
                        is the photocathode because it is the only one with
                        non-zero efficiency */
      {
        /* make sure the photon actually did hit the GEM stack */
        if (thePostPV->GetName() == "RICHHBDGEMFrontCu1Physical")
        {
          MakeHit(aStep);
        }
        break;
      }
      case FresnelReflection:
      case TotalInternalReflection:
      case LambertianReflection:
      case LobeReflection:
      case SpikeReflection:
      case BackScattering:
        fExpectedNextStatus = StepTooSmall;
        break;
      default:
        break;
      }
    }
  }

  return;
}

bool PHG4RICHSteppingAction::MakeHit(const G4Step* aStep)
{
  // collect energy and track
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  const G4Track* aTrack = aStep->GetTrack();
  const G4VTouchable* aTouch = aTrack->GetTouchable();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();

  // set sector number
  int sector_id = -999;
  bool sector_found = false;

  // Check if volume(1) is in sector volume, if not check volume(0)
  if (detector_->ePHENIXRICHConstruction::is_in_sector(aTouch->GetVolume(1)) > -1)
  {
    sector_id = aTouch->GetCopyNumber(1);
    sector_found = true;
  }
  else if (detector_->ePHENIXRICHConstruction::is_in_sector(aTouch->GetVolume()) > -1)
  {
    sector_id = aTouch->GetCopyNumber();
    sector_found = true;
  }

  if (!sector_found)
  {
    if (!aTouch->GetVolume(1) || !aTouch->GetVolume())
      cout << "WARNING: Missing volumes for hit!" << endl;
    else
      cout << "WARNING: Photon hit volume is not the RICH readout plane volume!" << endl;
  }

  hit = new PHG4Hitv1();

  //here we set the entrance values in cm
  hit->set_x(0, postPoint->GetPosition().x() / cm);
  hit->set_y(0, postPoint->GetPosition().y() / cm);
  hit->set_z(0, postPoint->GetPosition().z() / cm);

  // time in ns
  hit->set_t(0, postPoint->GetGlobalTime() / nanosecond);

  //same for exit values (photons absorbed/detected at boundary to post step volume)
  hit->set_x(1, postPoint->GetPosition().x() / cm);
  hit->set_y(1, postPoint->GetPosition().y() / cm);
  hit->set_z(1, postPoint->GetPosition().z() / cm);

  hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

  //set the track ID
  {
    hit->set_trkid(aTrack->GetTrackID());
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        hit->set_trkid(pp->GetUserTrackId());
        hit->set_shower_id(pp->GetShower()->get_id());

        pp->SetKeep(true);  // we want to keep the track
      }
    }
  }

  // set optical photon energy deposition
  hit->set_edep(edep);

  // Now add the hit
  hits_->AddHit(sector_id, hit);

  {
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        pp->GetShower()->add_g4hit_id(hits_->GetID(), hit->get_hit_id());
      }
    }
  }

  // return true to indicate the hit was used
  return true;
}

void PHG4RICHSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
