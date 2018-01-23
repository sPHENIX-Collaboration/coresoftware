#include "PHG4PSTOFSteppingAction.h"
#include "PHG4PSTOFDetector.h"
#include "PHG4Parameters.h"
#include "PHG4ParametersContainer.h"
#include "PHG4StepStatusDecode.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4PSTOFSteppingAction::PHG4PSTOFSteppingAction(PHG4PSTOFDetector* detector, const PHG4ParametersContainer* parameters)
  : detector_(detector)
  , paramscontainer(parameters)
  , hits_(nullptr)
  , hit(nullptr)
  , savehitcontainer(nullptr)
  , savevolpre(nullptr)
  , savevolpost(nullptr)
  , savetrackid(-1)
  , saveprestepstatus(-1)
  , savepoststepstatus(-1)
  , edepsum(0)
  , eionsum(0)
{}

PHG4PSTOFSteppingAction::~PHG4PSTOFSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4PSTOFSteppingAction::UserSteppingAction(const G4Step* aStep, bool was_used)
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
      aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
  {
    geantino = true;
  }
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  //       cout << "track id " << aTrack->GetTrackID() << endl;
  //       cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
  //       cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;

  layer_id = touch->GetCopyNumber();
  if (layer_id != whichactive)
  {
    cout << PHWHERE << " inconsistency between G4 copy number: " 
	 << layer_id << " and module id from detector: "
	 << whichactive << endl;
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
	cout << GetName() << ": New Hit for  " << endl;
cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
<< ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
           << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(saveprestepstatus)
           << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(savepoststepstatus) << endl;
      cout << "last track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID() << endl;
      cout << "phys pre vol: " << volume->GetName()
           << " post vol : " << touchpost->GetVolume()->GetName() << endl;
      cout << " previous phys pre vol: " << savevolpre->GetName()
           << " previous phys post vol: " << savevolpost->GetName() << endl;
      }
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
      cout << "implement stuff for whichactive < 0" << endl;
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
    if (!hit || !isfinite(hit->get_x(0)))
    {
      cout << GetName() << ": hit was not created" << endl;
      cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
           << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
           << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(saveprestepstatus)
           << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(savepoststepstatus) << endl;
      cout << "last track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID() << endl;
      cout << "phys pre vol: " << volume->GetName()
           << " post vol : " << touchpost->GetVolume()->GetName() << endl;
      cout << " previous phys pre vol: " << savevolpre->GetName()
           << " previous phys post vol: " << savevolpost->GetName() << endl;
      exit(1);
    }
    // check if track id matches the initial one when the hit was created
    if (aTrack->GetTrackID() != savetrackid)
    {
      cout << GetName() << ": hits do not belong to the same track" << endl;
      cout << "saved track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID()
           << ", prestep status: " << prePoint->GetStepStatus()
           << ", previous post step status: " << savepoststepstatus
           << endl;

      exit(1);
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
  string hitnodename;
  if (detector_->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + detector_->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + detector_->GetName();
  }

  //now look for the map and grab a pointer to it.
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());

  // if we do not find the node we need to make it.
  if (!hits_)
  {
    std::cout << "PHG4PSTOFSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
}
