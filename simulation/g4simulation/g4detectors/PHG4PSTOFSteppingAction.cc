#include "PHG4PSTOFSteppingAction.h"
#include "PHG4PSTOFDetector.h"
#include "PHG4StepStatusDecode.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>        // for G4TouchableHandle
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <cmath>  // for isfinite
#include <iostream>
#include <string>  // for operator<<, string

class PHCompositeNode;

//____________________________________________________________________________..
PHG4PSTOFSteppingAction::PHG4PSTOFSteppingAction(PHG4PSTOFDetector* detector, const PHParametersContainer* /*parameters*/)
  : PHG4SteppingAction(detector->GetName())
  , detector_(detector)
{
}

PHG4PSTOFSteppingAction::~PHG4PSTOFSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4PSTOFSteppingAction::UserSteppingAction(const G4Step* aStep, bool /*was_used*/)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();
  // IsInPSTOF(volume) returns
  //  == 0 outside of pstof
  //   > 0 for hits in active volume
  //  < 0 for hits in passive material
  int whichactive = detector_->IsInPSTOF(volume);
  if (!whichactive)
  {
    return false;
  }

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  const G4Track* aTrack = aStep->GetTrack();

  int layer_id = 0;  // what the heck is this?
  bool geantino = false;
  // the check for the pdg code speeds things up, I do not want to make
  // an expensive string compare for every track when we know
  // geantino or chargedgeantino has pid=0
  if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
      aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos)
  {
    geantino = true;
  }
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  //       std::cout << "track id " << aTrack->GetTrackID() << std::endl;
  //       std::cout << "time prepoint: " << prePoint->GetGlobalTime() << std::endl;
  //       std::cout << "time postpoint: " << postPoint->GetGlobalTime() << std::endl;

  layer_id = touch->GetCopyNumber();
  if (layer_id != whichactive)
  {
    std::cout << PHWHERE << " inconsistency between G4 copy number: "
              << layer_id << " and module id from detector: "
              << whichactive << std::endl;
    gSystem->Exit(1);
  }

  switch (prePoint->GetStepStatus())
  {
  case fPostStepDoItProc:
    if (savepoststepstatus != fGeomBoundary)
    {
      break;
    }
    else
    {
      std::cout << GetName() << ": New Hit for  " << std::endl;
      std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
                << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
                << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(saveprestepstatus)
                << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(savepoststepstatus) << std::endl;
      std::cout << "last track: " << savetrackid
                << ", current trackid: " << aTrack->GetTrackID() << std::endl;
      std::cout << "phys pre vol: " << volume->GetName()
                << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
      std::cout << " previous phys pre vol: " << savevolpre->GetName()
                << " previous phys post vol: " << savevolpost->GetName() << std::endl;
    }
    [[fallthrough]];
  case fGeomBoundary:
  case fUndefined:
    if (!hit)
    {
      hit = new PHG4Hitv1();
    }
    hit->set_layer(layer_id);
    //here we set the entrance values in cm
    hit->set_x(0, prePoint->GetPosition().x() / cm);
    hit->set_y(0, prePoint->GetPosition().y() / cm);
    hit->set_z(0, prePoint->GetPosition().z() / cm);
    // time in ns
    hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    //set the track ID
    hit->set_trkid(aTrack->GetTrackID());
    savetrackid = aTrack->GetTrackID();
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        hit->set_trkid(pp->GetUserTrackId());
      }
    }
    //set the initial energy deposit
    edepsum = 0;
    if (whichactive > 0)
    {
      eionsum = 0;
      hit->set_eion(0);
      savehitcontainer = hits_;
    }
    else
    {
      std::cout << "implement stuff for whichactive < 0" << std::endl;
      gSystem->Exit(1);
    }
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        hit->set_trkid(pp->GetUserTrackId());
        pp->GetShower()->add_g4hit_id(savehitcontainer->GetID(), hit->get_hit_id());
      }
    }

    break;
  default:
    break;
  }

  // some sanity checks for inconsistencies
  // check if this hit was created, if not print out last post step status
  if (!hit || !std::isfinite(hit->get_x(0)))
  {
    std::cout << GetName() << ": hit was not created" << std::endl;
    std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
              << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
              << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(saveprestepstatus)
              << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(savepoststepstatus) << std::endl;
    std::cout << "last track: " << savetrackid
              << ", current trackid: " << aTrack->GetTrackID() << std::endl;
    std::cout << "phys pre vol: " << volume->GetName()
              << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
    std::cout << " previous phys pre vol: " << savevolpre->GetName()
              << " previous phys post vol: " << savevolpost->GetName() << std::endl;
    gSystem->Exit(1);
  }
  // check if track id matches the initial one when the hit was created
  if (aTrack->GetTrackID() != savetrackid)
  {
    std::cout << GetName() << ": hits do not belong to the same track" << std::endl;
    std::cout << "saved track: " << savetrackid
              << ", current trackid: " << aTrack->GetTrackID()
              << ", prestep status: " << prePoint->GetStepStatus()
              << ", previous post step status: " << savepoststepstatus
              << std::endl;

    gSystem->Exit(1);
  }
  saveprestepstatus = prePoint->GetStepStatus();
  savepoststepstatus = postPoint->GetStepStatus();
  savevolpre = volume;
  savevolpost = touchpost->GetVolume();

  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  //sum up the energy to get total deposited
  edepsum += edep;
  if (whichactive > 0)
  {
    eionsum += eion;
  }
  // if any of these conditions is true this is the last step in
  // this volume and we need to save the hit
  // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
  // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
  // (happens when your detector goes outside world volume)
  // postPoint->GetStepStatus() == fAtRestDoItProc: track stops (typically
  // aTrack->GetTrackStatus() == fStopAndKill is also set)
  // aTrack->GetTrackStatus() == fStopAndKill: track ends
  if (postPoint->GetStepStatus() == fGeomBoundary ||
      postPoint->GetStepStatus() == fWorldBoundary ||
      postPoint->GetStepStatus() == fAtRestDoItProc ||
      aTrack->GetTrackStatus() == fStopAndKill)
  {
    // save only hits with energy deposit (or geantino)
    if (edepsum > 0 || geantino)
    {
      // update values at exit coordinates and set keep flag
      // of track to keep
      hit->set_x(1, postPoint->GetPosition().x() / cm);
      hit->set_y(1, postPoint->GetPosition().y() / cm);
      hit->set_z(1, postPoint->GetPosition().z() / cm);

      hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
      if (geantino)
      {
        hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
        if (whichactive > 0)
        {
          hit->set_eion(-1);
        }
      }
      else
      {
        hit->set_edep(edepsum);
      }
      if (whichactive > 0)
      {
        hit->set_eion(eionsum);
      }
      savehitcontainer->AddHit(layer_id, hit);
      // ownership has been transferred to container, set to null
      // so we will create a new hit for the next track
      hit = nullptr;
    }
    else
    {
      // if this hit has no energy deposit, just reset it for reuse
      // this means we have to delete it in the dtor. If this was
      // the last hit we processed the memory is still allocated
      hit->Reset();
    }
  }
  // return true to indicate the hit was used
  return true;
}

//____________________________________________________________________________..
void PHG4PSTOFSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  std::string hitnodename;
  if (detector_->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + detector_->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + detector_->GetName();
  }

  //now look for the map and grab a pointer to it.
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);

  // if we do not find the node we need to make it.
  if (!hits_)
  {
    std::cout << "PHG4PSTOFSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
}
