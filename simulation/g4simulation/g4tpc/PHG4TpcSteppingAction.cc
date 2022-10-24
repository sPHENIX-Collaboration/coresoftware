#include "PHG4TpcSteppingAction.h"
#include "PHG4TpcDetector.h"

#include <g4detectors/PHG4StepStatusDecode.h>

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fPostSt...
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

#include <cmath>    // for isfinite
#include <cstdlib>  // for exit
#include <iostream>
#include <string>  // for operator<<, operator+

class PHCompositeNode;

//____________________________________________________________________________..
PHG4TpcSteppingAction::PHG4TpcSteppingAction(PHG4TpcDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_IsBlackHoleFlag(m_Params->get_int_param("blackhole"))
{
  if (std::isfinite(m_Params->get_double_param("steplimits")))
  {
    m_UseG4StepsFlag = 1;
  }
  SetName(m_Detector->GetName());
}

PHG4TpcSteppingAction::~PHG4TpcSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}
//____________________________________________________________________________..
bool PHG4TpcSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  // m_Detector->IsInTpc(volume)
  // returns
  //  0 is outside of Tpc
  //  1 is inside tpc gas
  // <0 is in tpc support structure (cage, endcaps,...)

  int whichactive = m_Detector->IsInTpc(volume);
  if (!whichactive)
  {
    return false;
  }
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (m_IsBlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

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
  int prepointstatus = prePoint->GetStepStatus();

  //       std::cout << "track id " << aTrack->GetTrackID() << std::endl;
  //       std::cout << "time prepoint: " << prePoint->GetGlobalTime() << std::endl;
  //       std::cout << "time postpoint: " << postPoint->GetGlobalTime() << std::endl;

  if ((m_UseG4StepsFlag > 0 && whichactive > 0) ||
      prepointstatus == fGeomBoundary ||
      prepointstatus == fUndefined ||
      (prepointstatus == fPostStepDoItProc && m_SavePostStepStatus == fGeomBoundary))
  {
    unsigned int layer_id = 99;  // no layer number for the hit, use a non-existent one for now, replace it later
    // this is for debugging weird occurances we have occasionally
    if (prepointstatus == fPostStepDoItProc && m_SavePostStepStatus == fGeomBoundary)
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
    // if previous hit was saved, hit pointer was set to nullptr
    // and we have to make a new one
    if (!m_Hit)
    {
      m_Hit = new PHG4Hitv1();
    }
    m_Hit->set_layer(layer_id);
    //here we set the entrance values in cm
    m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
    m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
    m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

    // momentum
    m_Hit->set_px(0, prePoint->GetMomentum().x() / GeV);
    m_Hit->set_py(0, prePoint->GetMomentum().y() / GeV);
    m_Hit->set_pz(0, prePoint->GetMomentum().z() / GeV);

    // time in ns
    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    //set and save the track ID
    m_Hit->set_trkid(aTrack->GetTrackID());
    m_SaveTrackId = aTrack->GetTrackID();
    //set the initial energy deposit
    m_Hit->set_edep(0);
    if (whichactive > 0)  // return of IsInTpcDetector, > 0 hit in tpc gas volume, < 0 hit in support structures
    {
      m_Hit->set_eion(0);
      // Now save the container we want to add this hit to
      m_CurrentHitContainer = m_HitContainer;
    }
    else
    {
      m_CurrentHitContainer = m_AbsorberHitContainer;
    }
    if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        m_Hit->set_trkid(pp->GetUserTrackId());
        m_Hit->set_shower_id(pp->GetShower()->get_id());
        m_Shower = pp->GetShower();
      }
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
      exit(1);
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

      exit(1);
    }
    m_SavePreStepStatus = prePoint->GetStepStatus();
    m_SavePostStepStatus = postPoint->GetStepStatus();
    m_SaveVolPre = volume;
    m_SaveVolPost = touchpost->GetVolume();
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    m_Hit->set_px(1, postPoint->GetMomentum().x() / GeV);
    m_Hit->set_py(1, postPoint->GetMomentum().y() / GeV);
    m_Hit->set_pz(1, postPoint->GetMomentum().z() / GeV);

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    //sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (whichactive > 0)
    {
      m_Hit->set_eion(m_Hit->get_eion() + eion);
    }
    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      if (whichactive > 0)
      {
        m_Hit->set_eion(-1);
      }
    }
    if (edep > 0)
    {
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          pp->SetKeep(1);  // we want to keep the track
        }
      }
    }

    // if any of these conditions is true this is the last step in
    // this volume and we need to save the hit
    // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
    // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
    // (happens when your detector goes outside world volume)
    // postPoint->GetStepStatus() == fAtRestDoItProc: track stops (typically
    // aTrack->GetTrackStatus() == fStopAndKill is also set)
    // aTrack->GetTrackStatus() == fStopAndKill: track ends
    if ((m_UseG4StepsFlag > 0 && whichactive > 0) ||
        postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill)
    {
      // save only hits with energy deposit (or -1 for geantino)
      if (m_Hit->get_edep())
      {
        m_CurrentHitContainer->AddHit(layer_id, m_Hit);
        if (m_Shower)
        {
          m_Shower->add_g4hit_id(m_CurrentHitContainer->GetID(), m_Hit->get_hit_id());
        }

        double rin = sqrt(m_Hit->get_x(0) * m_Hit->get_x(0) + m_Hit->get_y(0) * m_Hit->get_y(0));
        double rout = sqrt(m_Hit->get_x(1) * m_Hit->get_x(1) + m_Hit->get_y(1) * m_Hit->get_y(1));
        if (Verbosity() > 10)
          if ((rin > 69.0 && rin < 70.125) || (rout > 69.0 && rout < 70.125))
          {
            std::cout << "Added Tpc g4hit with rin, rout = " << rin << "  " << rout
                      << " g4hitid " << m_Hit->get_hit_id() << std::endl;
            std::cout << " xin " << m_Hit->get_x(0)
                      << " yin " << m_Hit->get_y(0)
                      << " zin " << m_Hit->get_z(0)
                      << " rin " << rin
                      << std::endl;
            std::cout << " xout " << m_Hit->get_x(1)
                      << " yout " << m_Hit->get_y(1)
                      << " zout " << m_Hit->get_z(1)
                      << " rout " << rout
                      << std::endl;
            std::cout << " xav " << (m_Hit->get_x(1) + m_Hit->get_x(0)) / 2.0
                      << " yav " << (m_Hit->get_y(1) + m_Hit->get_y(0)) / 2.0
                      << " zav " << (m_Hit->get_z(1) + m_Hit->get_z(0)) / 2.0
                      << " rav " << (rout + rin) / 2.0
                      << std::endl;
          }

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
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4TpcSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_AbsorberNodeName);

  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4TpcSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
  }
  if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 1)
    {
      std::cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << m_AbsorberNodeName << std::endl;
    }
  }
}

void PHG4TpcSteppingAction::SetHitNodeName(const std::string& type, const std::string& name)
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
