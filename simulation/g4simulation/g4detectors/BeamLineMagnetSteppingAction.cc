#include "BeamLineMagnetSteppingAction.h"

#include "BeamLineMagnetDetector.h"
#include "PHG4StepStatusDecode.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

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
BeamLineMagnetSteppingAction::BeamLineMagnetSteppingAction(BeamLineMagnetDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_AbsorberActiveFlag(m_Params->get_int_param("absorberactive"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
{
}

BeamLineMagnetSteppingAction::~BeamLineMagnetSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool BeamLineMagnetSteppingAction::UserSteppingAction(const G4Step* aStep, bool /*was_used*/)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();
  // IsInVTX(volume) returns
  //  == 0 outside of pstof
  //   > 0 for hits in active volume
  //  < 0 for hits in passive material
  int whichactive = m_Detector->IsInBeamLineMagnet(volume);
  if (!whichactive)
  {
    return false;
  }

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything ONLY if in the magnet iron, just put all kinetic energy into edep
  if (whichactive == -1)
  {
    if (m_BlackHoleFlag)
    {
      edep = aTrack->GetKineticEnergy() / GeV;
      G4Track* killtrack = const_cast<G4Track*>(aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }
    if (!m_AbsorberActiveFlag)
    {
      return false;
    }
  }
  if (whichactive == -2)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }
  if (!m_ActiveFlag)
  {
    return false;
  }

  int magnet_id = volume->GetCopyNo();  // magnet id is stored in copy number
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

  switch (prePoint->GetStepStatus())
  {
  case fPostStepDoItProc:
    if (m_SavePostStepStatus != fGeomBoundary)
    {
      break;
    }
    else
    {
      std::cout << GetName() << ": New Hit for  " << std::endl;
      std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
                << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
                << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
                << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
      std::cout << "last track: " << m_SaveTrackId
                << ", current trackid: " << aTrack->GetTrackID() << std::endl;
      std::cout << "phys pre vol: " << volume->GetName()
                << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
      std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
                << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    }
    [[fallthrough]];
  case fGeomBoundary:
  case fUndefined:
    if (!m_Hit)
    {
      m_Hit = new PHG4Hitv1();
    }
    m_Hit->set_layer(magnet_id);
    //here we set the entrance values in cm
    m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
    m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
    m_Hit->set_z(0, prePoint->GetPosition().z() / cm);
    // time in ns
    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    //set the track ID
    m_Hit->set_trkid(aTrack->GetTrackID());
    m_SaveTrackId = aTrack->GetTrackID();
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
      }
    }
    //set the initial energy deposit
    m_EdepSum = 0;
    if (whichactive > 0)
    {
      m_SaveHitContainer = m_HitContainer;
    }
    else
    {
      m_SaveHitContainer = m_AbsorberHitContainer;
    }
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
        pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());
      }
    }

    break;
  default:
    break;
  }

  // some sanity checks for inconsistencies
  // check if this hit was created, if not print out last post step status
  if (!m_Hit || !std::isfinite(m_Hit->get_x(0)))
  {
    std::cout << GetName() << ": hit was not created" << std::endl;
    std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
              << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
              << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
              << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
    std::cout << "last track: " << m_SaveTrackId
              << ", current trackid: " << aTrack->GetTrackID() << std::endl;
    std::cout << "phys pre vol: " << volume->GetName()
              << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
    std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
              << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    gSystem->Exit(1);
  }
  // check if track id matches the initial one when the hit was created
  if (aTrack->GetTrackID() != m_SaveTrackId)
  {
    std::cout << GetName() << ": hits do not belong to the same track" << std::endl;
    std::cout << "saved track: " << m_SaveTrackId
              << ", current trackid: " << aTrack->GetTrackID()
              << ", prestep status: " << prePoint->GetStepStatus()
              << ", previous post step status: " << m_SavePostStepStatus
              << std::endl;

    gSystem->Exit(1);
  }
  m_SavePreStepStatus = prePoint->GetStepStatus();
  m_SavePostStepStatus = postPoint->GetStepStatus();
  m_SaveVolPre = volume;
  m_SaveVolPost = touchpost->GetVolume();

  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  //sum up the energy to get total deposited
  m_EdepSum += edep;
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
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
      if (geantino)
      {
        m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      }
      else
      {
        m_Hit->set_edep(m_EdepSum);
      }
      m_SaveHitContainer->AddHit(magnet_id, m_Hit);
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
void BeamLineMagnetSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_AbsorberNodeName);

  // If this is a black hole but otherwise not an active detector there is
  // no hitnode, otherwise this will crash in the user stepping action
  if (!m_HitContainer)
  {
    if (!m_BlackHoleFlag)
    {
      std::cout << "BeamLineMagnetSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
    }
  }
  if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 1)
    {
      std::cout << "BeamLineMagnetSteppingAction::SetTopNode - unable to find " << m_AbsorberNodeName << std::endl;
    }
  }
}

void BeamLineMagnetSteppingAction::SetHitNodeName(const std::string& type, const std::string& name)
{
  if (type == "G4HIT")
  {
    m_HitNodeName = name;
    return;
  }
  else if (type == "G4HIT_ABSORBER")
  {
    m_AbsorberNodeName = name;
    return;
  }
  std::cout << "Invalid output hit node type " << type << std::endl;
  gSystem->Exit(1);
  return;
}
