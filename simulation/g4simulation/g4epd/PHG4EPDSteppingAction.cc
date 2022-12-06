/* vim: set sw=2 ft=cpp: */

#include "PHG4EPDSteppingAction.h"

#include "PHG4EPDDetector.h"

#include <g4detectors/PHG4StepStatusDecode.h>

#include <phool/getClass.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <TSystem.h>

#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ReferenceCountedHandle.hh>
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>
#include <Geant4/G4StepStatus.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackStatus.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4VTouchable.hh>
#include <Geant4/G4VUserTrackInformation.hh>

#include <cstdint>  // for int32_t
#include <iostream>
#include <string>

class G4VPhysicalVolume;

PHG4EPDSteppingAction::PHG4EPDSteppingAction(PHG4EPDDetector* detector,
                                             const PHParameters* /*unused*/)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
{
}

PHG4EPDSteppingAction::~PHG4EPDSteppingAction()
{
  delete m_Hit;
}

bool PHG4EPDSteppingAction::UserSteppingAction(const G4Step* aStep, bool /*was_used*/)
{
  G4StepPoint* prestep = aStep->GetPreStepPoint();
  G4TouchableHandle prehandle = prestep->GetTouchableHandle();

  G4VPhysicalVolume* volume = prehandle->GetVolume();

  // m_Detector->IsInDetector(volume)
  // returns
  //  0 is outside of EPD
  //  1 is inside scintillator
  // -1 is inside support structure (dead material)
  int whichactive = m_Detector->IsInDetector(volume);
  if (!whichactive)
  {
    return false;
  }

  G4double deposit = aStep->GetTotalEnergyDeposit() / GeV;
  G4double ionising = deposit - aStep->GetNonIonizingEnergyDeposit() / GeV;
  G4double light_yield = GetVisibleEnergyDeposition(aStep) / GeV;

  //  G4StepStatus prestatus = prestep->GetStepStatus();

  int32_t tile_id = m_Detector->module_id_for(volume);

  G4Track const* track = aStep->GetTrack();

  G4ParticleDefinition const* particle = track->GetParticleDefinition();

  bool geantino = (particle->GetPDGEncoding() == 0 && particle->GetParticleName().find("geantino") != std::string::npos);

  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();

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
                << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << std::endl;
    }
    [[fallthrough]];
  case fGeomBoundary:
  case fUndefined:
    if (m_Hit == nullptr)
    {
      m_Hit = new PHG4Hitv1();
    }

    // only for active columes (scintillators)
    if (whichactive > 0)
    {
      m_Hit->set_eion(0);
      m_Hit->set_light_yield(0);
    }
    m_Hit->set_x(0, prestep->GetPosition().x() / cm);
    m_Hit->set_y(0, prestep->GetPosition().y() / cm);
    m_Hit->set_z(0, prestep->GetPosition().z() / cm);
    m_Hit->set_t(0, prestep->GetGlobalTime() / nanosecond);

    m_Hit->set_trkid(track->GetTrackID());

    if (PHG4TrackUserInfoV1* userinfo = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()))
    {
      m_Hit->set_trkid(userinfo->GetUserTrackId());

      userinfo->GetShower()->add_g4hit_id(m_HitContainer->GetID(), m_Hit->get_hit_id());
    }

    m_Hit->set_edep(0);
    break;
  default:
    break;
  }

  G4StepPoint* poststep = aStep->GetPostStepPoint();
  const G4ThreeVector postpos = poststep->GetPosition();
  m_SavePostStepStatus = postPoint->GetStepStatus();
  m_Hit->set_edep(m_Hit->get_edep() + deposit);
  if (whichactive > 0)
  {
    m_Hit->set_eion(m_Hit->get_eion() + ionising);
    m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield * GetLightCorrection(postpos.x(), postpos.y()));
  }

  if (postPoint->GetStepStatus() != fGeomBoundary &&
      postPoint->GetStepStatus() != fWorldBoundary &&
      postPoint->GetStepStatus() != fAtRestDoItProc &&
      track->GetTrackStatus() != fStopAndKill)
  {
    return true;
  }
  if (m_Hit->get_edep() <= 0 && !geantino)
  {
    m_Hit->Reset();

    return true;
  }

  m_Hit->set_x(1, poststep->GetPosition().x() / cm);
  m_Hit->set_y(1, poststep->GetPosition().y() / cm);
  m_Hit->set_z(1, poststep->GetPosition().z() / cm);
  m_Hit->set_t(1, poststep->GetGlobalTime() / nanosecond);

  PHG4TrackUserInfoV1* userinfo = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation());

  if (userinfo != nullptr)
  {
    userinfo->SetKeep(1);
  }

  if (geantino)
  {
    m_Hit->set_edep(-1.);
    m_Hit->set_eion(-1.);
    m_Hit->set_light_yield(-1.);
  }

  m_HitContainer->AddHit(tile_id, m_Hit);

  m_Hit = nullptr;

  return true;
}

void PHG4EPDSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);

  m_SupportHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_SupportNodeName);
  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4EPDSteppingAction::SetTopNode - unable to find hit node " << m_HitNodeName << std::endl;
    gSystem->Exit(1);
  }
  // this is perfectly fine if support hits are disabled
  if (!m_SupportHitContainer)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4EPDSteppingAction::SetTopNode - unable to find support hit node " << m_SupportNodeName << std::endl;
    }
  }
}

void PHG4EPDSteppingAction::SetHitNodeName(const std::string& type, const std::string& name)
{
  if (type == "G4HIT")
  {
    m_HitNodeName = name;
    return;
  }
  else if (type == "G4HIT_SUPPORT")
  {
    m_SupportNodeName = name;
    return;
  }
  std::cout << "Invalid output hit node type " << type << std::endl;
  gSystem->Exit(1);
  return;
}
