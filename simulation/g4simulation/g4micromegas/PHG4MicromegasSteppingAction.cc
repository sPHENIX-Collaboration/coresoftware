/*!
 * \file PHG4MicromegasSteppingAction.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasSteppingAction.h"

#include "PHG4MicromegasDetector.h"

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
PHG4MicromegasSteppingAction::PHG4MicromegasSteppingAction(PHG4MicromegasDetector *detector, const PHParameters *parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
{}

//____________________________________________________________________________..
// This is the implementation of the G4 UserSteppingAction
bool PHG4MicromegasSteppingAction::UserSteppingAction(const G4Step *aStep,bool /*was_used*/)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();

  // get volume of the current step
  G4VPhysicalVolume *volume = touch->GetVolume();
  if (!m_Detector->IsInDetector(volume)) return false;

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  const G4Track *aTrack = aStep->GetTrack();

  // if this detector stops everything, just put all kinetic energy into edep
  if (m_BlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    auto killtrack = const_cast<G4Track *>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  bool geantino =
    aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
    aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos;

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
  switch (prePoint->GetStepStatus())
  {

    case fPostStepDoItProc:
    if (m_SavePostStepStatus == fGeomBoundary)
    {

      // this is an impossible G4 Step print out diagnostic to help debug, not sure if
      // this is still with us
      std::cout << "PHG4MicromegasSteppingAction::UserSteppingAction - " << GetName() << ": New Hit for  " << std::endl;
      std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
        << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
        << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
        << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus)
        << std::endl;

      std::cout << "last track: " << m_SaveTrackId << ", current trackid: " << aTrack->GetTrackID() << std::endl;
      std::cout << "phys pre vol: " << volume->GetName() << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
      std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName() << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;
    }
    break;

    // These are the normal cases
    case fGeomBoundary:
    case fUndefined:
    {
      if (!m_hit) m_hit.reset( new PHG4Hitv1() );

      // assign layer
      m_hit->set_layer(m_Detector->get_layer(volume));

      const auto tileid = m_Detector->get_tileid(volume);
      m_hit->set_property(PHG4Hit::prop_index_i, tileid);

      // here we set the entrance values in cm
      m_hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_hit->set_z(0, prePoint->GetPosition().z() / cm);

      // momentum
      m_hit->set_px(0, prePoint->GetMomentum().x() / GeV);
      m_hit->set_py(0, prePoint->GetMomentum().y() / GeV);
      m_hit->set_pz(0, prePoint->GetMomentum().z() / GeV);

      // time in ns
      m_hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);

      // set the track ID
      m_hit->set_trkid(aTrack->GetTrackID());
      m_SaveTrackId = aTrack->GetTrackID();

      // reset the initial energy deposit
      m_EdepSum = 0;
      m_EionSum = 0;
      m_hit->set_edep(0);
      m_hit->set_eion(0);
      m_SaveHitContainer = m_hitContainer;

      // this is for the tracking of the truth info
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
        {
          m_hit->set_trkid(pp->GetUserTrackId());
          pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(), m_hit->get_hit_id());
        }
      }
      break;
    }

    default: break;

  }

  if (!m_hit || !std::isfinite(m_hit->get_x(0)))
  {
    std::cout << GetName() << ": hit was not created" << std::endl;
    std::cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
      << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
      << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
      << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus)
      << std::endl;
    std::cout << "last track: " << m_SaveTrackId << ", current trackid: " << aTrack->GetTrackID() << std::endl;
    std::cout << "phys pre vol: " << volume->GetName() << " post vol : " << touchpost->GetVolume()->GetName() << std::endl;
    std::cout << " previous phys pre vol: " << m_SaveVolPre->GetName() << " previous phys post vol: " << m_SaveVolPost->GetName() << std::endl;

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
  m_EionSum += eion;

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
      m_hit->set_x(1, postPoint->GetPosition().x() / cm);
      m_hit->set_y(1, postPoint->GetPosition().y() / cm);
      m_hit->set_z(1, postPoint->GetPosition().z() / cm);

      m_hit->set_px(1, postPoint->GetMomentum().x() / GeV);
      m_hit->set_py(1, postPoint->GetMomentum().y() / GeV);
      m_hit->set_pz(1, postPoint->GetMomentum().z() / GeV);

      m_hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
      if (G4VUserTrackInformation *p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1 *pp = dynamic_cast<PHG4TrackUserInfoV1 *>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
      if (geantino)
      {
        m_hit->set_edep(-1);
        m_hit->set_eion(-1);
      }
      else
      {
        m_hit->set_edep(m_EdepSum);
        m_hit->set_eion(m_EionSum);
      }


      // add in container
      m_SaveHitContainer->AddHit(m_hit->get_layer(), m_hit.get());

      // ownership is transferred to container
      // so we release the hit
      m_hit.release();

    } else {

      m_hit->Reset();

    }
  }
  // return true to indicate the hit was used
  return true;
}

//____________________________________________________________________________..
void PHG4MicromegasSteppingAction::SetInterfacePointers(PHCompositeNode *topNode)
{
  std::string hitnodename = "G4HIT_" + m_Detector->SuperDetector();
  m_hitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!m_hitContainer) std::cout << "PHG4MicromegasSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
}
