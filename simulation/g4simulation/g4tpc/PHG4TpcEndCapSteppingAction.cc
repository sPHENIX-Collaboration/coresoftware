#include "PHG4TpcEndCapSteppingAction.h"

#include "PHG4TpcEndCapDetector.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4StepStatusDecode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ReferenceCountedHandle.hh>
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>
#include <Geant4/G4StepStatus.hh>
#include <Geant4/G4String.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackStatus.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4VTouchable.hh>
#include <Geant4/G4VUserTrackInformation.hh>

#include <cmath>
#include <iostream>
#include <string>

class PHCompositeNode;

//____________________________________________________________________________..
PHG4TpcEndCapSteppingAction::PHG4TpcEndCapSteppingAction(PHG4TpcEndCapDetector *detector, const PHParameters *parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
{
}

//____________________________________________________________________________..
PHG4TpcEndCapSteppingAction::~PHG4TpcEndCapSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
// This is the implementation of the G4 UserSteppingAction
bool PHG4TpcEndCapSteppingAction::UserSteppingAction(const G4Step *aStep, bool /*was_used*/)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume *volume = touch->GetVolume();
  // IsInDetector(volume) returns
  //  == 0 outside of detector
  //   > 0 for hits in active volume
  //  < 0 for hits in passive material
  int whichactive = m_Detector->IsInDetector(volume);
  if (!whichactive)
  {
    return false;
  }
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  const G4Track *aTrack = aStep->GetTrack();
  // if this detector stops everything, just put all kinetic energy into edep
  if (m_BlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track *killtrack = const_cast<G4Track *>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }
  // we use here only one detector in this simple example
  // if you deal with multiple detectors in this stepping action
  // the detector id can be used to distinguish between them
  // hits can easily be analyzed later according to their detector id
  int detector_id = 0;  // we use here only one detector in this simple example
  bool geantino = false;
  // the check for the pdg code speeds things up, I do not want to make
  // an expensive string compare for every track when we know
  // geantino or chargedgeantino has pid=0
  if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
      aTrack->GetParticleDefinition()->GetParticleName().find("geantino") !=
          std::string::npos)  // this also accounts for "chargedgeantino"
  {
    geantino = true;
  }
  G4StepPoint *prePoint = aStep->GetPreStepPoint();
  G4StepPoint *postPoint = aStep->GetPostStepPoint();

  // Here we have to decide if we need to create a new hit.  Normally this should
  // only be neccessary if a G4 Track enters a new volume or is freshly created
  // For this we look at the step status of the prePoint (beginning of the G4 Step).
  // This should be either fGeomBoundary (G4 Track crosses into volume) or
  // fUndefined (G4 Track newly created)
  // Sadly over the years with different G4 versions we have observed cases where
  // G4 produces "impossible hits" which we try to catch here
  // These errors were always rare and it is not clear if they still exist but we
  // still check for them for safety. We can reproduce G4 runs identically (if given
  // the sequence of random number seeds you find in the log), the printouts help
  // us giving the G4 support information about those failures
  //
  switch (prePoint->GetStepStatus())
  {
  case fPostStepDoItProc:
    if (m_SavePostStepStatus != fGeomBoundary)
    {
      // this is the okay case, fPostStepDoItProc called in a volume, not first thing inside
      // a new volume, just proceed here
      break;
    }
    else
    {
      // this is an impossible G4 Step print out diagnostic to help debug, not sure if
      // this is still with us
      std::cout << GetName() << ": New Hit for  " << std::endl;
      std::cout << "prestep status: "
                << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
                << ", poststep status: "
                << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
                << ", last pre step status: "
                << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
                << ", last post step status: "
                << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
      std::cout << "last track: " << m_SaveTrackId
                << ", current trackid: " << aTrack->GetTrackID() << std::endl;
      std::cout << "phys pre vol: " << volume->GetName()
                << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
      std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
                << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    }
    [[fallthrough]];
    // These are the normal cases
  case fGeomBoundary:
  case fUndefined:
    if (!m_Hit)
    {
      m_Hit = new PHG4Hitv1();
    }
    m_Hit->set_layer(detector_id);
    // here we set the entrance values in cm
    m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
    m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
    m_Hit->set_z(0, prePoint->GetPosition().z() / cm);
    // time in ns
    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    // set the track ID
    m_Hit->set_trkid(aTrack->GetTrackID());
    m_SaveTrackId = aTrack->GetTrackID();
    // set the initial energy deposit
    m_EdepSum = 0;
    // implement your own here://
    // add the properties you are interested in via set_XXX methods
    // you can find existing set methods in $OFFLINE_MAIN/include/g4main/PHG4Hit.h
    // this is initialization of your value. This is not needed you can just set the final
    // value at the last step in this volume later one
    if (whichactive > 0)
    {
      m_EionSum = 0;  // assuming the ionization energy is only needed for active
                      // volumes (scintillators)
      m_Hit->set_eion(0);
      m_SaveHitContainer = m_HitContainer;
    }
    else
    {
      std::cout << "implement stuff for whichactive < 0 (inactive volumes)" << std::endl;
      gSystem->Exit(1);
    }
    // this is for the tracking of the truth info
    if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
        pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());
      }
    }
    break;
  default:
    break;
  }

  // This section is called for every step
  // some sanity checks for inconsistencies (aka bugs) we have seen over the years
  // check if this hit was created, if not print out last post step status
  if (!m_Hit || !std::isfinite(m_Hit->get_x(0)))
  {
    std::cout << GetName() << ": hit was not created" << std::endl;
    std::cout << "prestep status: "
              << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
              << ", poststep status: "
              << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
              << ", last pre step status: "
              << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
              << ", last post step status: "
              << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
    std::cout << "last track: " << m_SaveTrackId
              << ", current trackid: " << aTrack->GetTrackID() << std::endl;
    std::cout << "phys pre vol: " << volume->GetName()
              << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
    std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
              << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    // This is fatal - a hit from nowhere. This needs to be looked at and fixed
    gSystem->Exit(1);
  }
  // check if track id matches the initial one when the hit was created
  if (aTrack->GetTrackID() != m_SaveTrackId)
  {
    std::cout << GetName() << ": hits do not belong to the same track" << std::endl;
    std::cout << "saved track: " << m_SaveTrackId
              << ", current trackid: " << aTrack->GetTrackID()
              << ", prestep status: " << prePoint->GetStepStatus()
              << ", previous post step status: " << m_SavePostStepStatus << std::endl;
    // This is fatal - a hit from nowhere. This needs to be looked at and fixed
    gSystem->Exit(1);
  }

  // We need to cache a few things from one step to the next
  // to identify impossible hits and subsequent debugging printout
  m_SavePreStepStatus = prePoint->GetStepStatus();
  m_SavePostStepStatus = postPoint->GetStepStatus();
  m_SaveVolPre = volume;
  m_SaveVolPost = touchpost->GetVolume();
  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  // sum up the energy to get total deposited
  m_EdepSum += edep;
  if (whichactive > 0)
  {
    m_EionSum += eion;
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
    if (m_EdepSum > 0 || geantino)
    {
      // update values at exit coordinates and set keep flag
      // of track to keep
      m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
      m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
      m_Hit->set_z(1, postPoint->GetPosition().z() / cm);
      m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
      if (geantino)
      {
        //implement your own here://
        // if you want to do something special for geantinos (normally you do not)
        m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way
                              // geantinos survive the g4hit compression
        if (whichactive > 0)
        {
          m_Hit->set_eion(-1);
        }
      }
      else
      {
        m_Hit->set_edep(m_EdepSum);
      }
      //implement your own here://
      // what you set here will be saved in the output
      if (whichactive > 0)
      {
        m_Hit->set_eion(m_EionSum);
      }
      m_SaveHitContainer->AddHit(detector_id, m_Hit);
      // ownership has been transferred to container, set to null
      // so we will create a new hit for the next track
      m_Hit = nullptr;
    }
    else
    {
      // if this hit has no energy deposit, just reset it for reuse
      // this means we have to delete it in the dtor. If this was
      // the last hit we processed the memory is still allocated
      m_Hit->Reset();
    }
  }
  // return true to indicate the hit was used
  return true;
}

//____________________________________________________________________________..
void PHG4TpcEndCapSteppingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  // now look for the map and grab a pointer to it.
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  // if we do not find the node we need to make it.
  if (!m_HitContainer)
  {
    std::cout << "PHG4TpcEndCapSteppingAction::SetTopNode - unable to find "
              << m_HitNodeName << std::endl;
    gSystem->Exit(1);
  }
}

void PHG4TpcEndCapSteppingAction::SetHitNodeName(const std::string &type, const std::string &name)
{
  if (type == "G4HIT")
  {
    m_HitNodeName = name;
    return;
  }
  std::cout << "Invalid output hit node type " << type << std::endl;
  gSystem->Exit(1);
  return;
}
