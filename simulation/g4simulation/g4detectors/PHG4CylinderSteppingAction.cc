#include "PHG4CylinderSteppingAction.h"
#include "PHG4CylinderDetector.h"
#include "PHG4CylinderSubsystem.h"

#include "PHG4StepStatusDecode.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/io/ios_state.hpp>
#pragma GCC diagnostic pop

#include <cmath>    // for isfinite, copysign
#include <cstdlib>  // for exit
#include <iomanip>
#include <iostream>
#include <string>  // for operator<<, char_traits

class PHCompositeNode;

//____________________________________________________________________________..
PHG4CylinderSteppingAction::PHG4CylinderSteppingAction(PHG4CylinderSubsystem* subsys, PHG4CylinderDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Subsystem(subsys)
  , m_Detector(detector)
  , m_Params(parameters)
  , m_HitContainer(nullptr)
  , m_Hit(nullptr)
  , m_SaveShower(nullptr)
  , m_SaveVolPre(nullptr)
  , m_SaveVolPost(nullptr)
  , m_SaveLightYieldFlag(m_Params->get_int_param("lightyield"))
  , m_SaveTrackId(-1)
  , m_SavePreStepStatus(-1)
  , m_SavePostStepStatus(-1)
  , m_ActiveFlag(m_Params->get_int_param("active"))
  , m_BlackHoleFlag(m_Params->get_int_param("blackhole"))
  , m_UseG4StepsFlag(m_Params->get_int_param("use_g4steps"))
  , m_Zmin(m_Params->get_double_param("place_z") * cm - m_Params->get_double_param("length") * cm / 2.)
  , m_Zmax(m_Params->get_double_param("place_z") * cm + m_Params->get_double_param("length") * cm / 2.)
  , m_Tmin(m_Params->get_double_param("tmin") * ns)
  , m_Tmax(m_Params->get_double_param("tmax") * ns)
{
  // G4 seems to have issues in the um range
  m_Zmin -= copysign(m_Zmin, 1. / 1e6 * cm);
  m_Zmax += copysign(m_Zmax, 1. / 1e6 * cm);
}

PHG4CylinderSteppingAction::~PHG4CylinderSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4CylinderSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  // get volume of the current step
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();

  G4VPhysicalVolume* volume = touch->GetVolume();
  // G4 just calls  UserSteppingAction for every step (and we loop then over all our
  // steppingactions. First we have to check if we are actually in our volume
  if (!m_Detector->IsInCylinder(volume))
  {
    return false;
  }

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;

  const G4Track* aTrack = aStep->GetTrack();

  // if this cylinder stops everything, just put all kinetic energy into edep
  if (m_BlackHoleFlag)
  {
    if ((!std::isfinite(m_Tmin) && !std::isfinite(m_Tmax)) ||
        aTrack->GetGlobalTime() < m_Tmin ||
        aTrack->GetGlobalTime() > m_Tmax)
    {
      edep = aTrack->GetKineticEnergy() / GeV;
      G4Track* killtrack = const_cast<G4Track*>(aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }
  }

  int layer_id = m_Detector->get_Layer();
  // test if we are active
  if (m_ActiveFlag)
  {
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
    //        std::cout << "time prepoint: " << prePoint->GetGlobalTime()/ns << std::endl;
    //        std::cout << "time postpoint: " << postPoint->GetGlobalTime()/ns << std::endl;
    //        std::cout << "kinetic energy: " <<  aTrack->GetKineticEnergy()/GeV << std::endl;
    //       G4ParticleDefinition* def = aTrack->GetDefinition();
    //       std::cout << "Particle: " << def->GetParticleName() << std::endl;
    int prepointstatus = prePoint->GetStepStatus();
    if (prepointstatus == fGeomBoundary ||
        prepointstatus == fUndefined ||
        (prepointstatus == fPostStepDoItProc && m_SavePostStepStatus == fGeomBoundary) ||
        m_UseG4StepsFlag > 0)
    {
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

      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }

      m_Hit->set_layer((unsigned int) layer_id);

      //here we set the entrance values in cm
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

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
      if (!geantino && !m_BlackHoleFlag)
      {
        m_Hit->set_eion(0);
      }
      if (m_SaveLightYieldFlag)
      {
        m_Hit->set_light_yield(0);
      }
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          m_Hit->set_trkid(pp->GetUserTrackId());
          m_Hit->set_shower_id(pp->GetShower()->get_id());
          m_SaveShower = pp->GetShower();
        }
      }

      if (!hasMotherSubsystem() && (m_Hit->get_z(0) * cm > m_Zmax || m_Hit->get_z(0) * cm < m_Zmin))
      {
        boost::io::ios_precision_saver ips(std::cout);
        std::cout << m_Detector->SuperDetector() << std::setprecision(9)
                  << " PHG4CylinderSteppingAction: Entry hit z " << m_Hit->get_z(0) * cm
                  << " outside acceptance,  zmin " << m_Zmin
                  << ", zmax " << m_Zmax << ", layer: " << layer_id << std::endl;
      }
    }
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
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
    m_SavePostStepStatus = postPoint->GetStepStatus();
    // check if track id matches the initial one when the hit was created
    if (aTrack->GetTrackID() != m_SaveTrackId)
    {
      std::cout << "hits do not belong to the same track" << std::endl;
      std::cout << "saved track: " << m_SaveTrackId
                << ", current trackid: " << aTrack->GetTrackID()
                << std::endl;
      exit(1);
    }
    m_SavePreStepStatus = prePoint->GetStepStatus();
    m_SavePostStepStatus = postPoint->GetStepStatus();
    m_SaveVolPre = volume;
    m_SaveVolPost = touchpost->GetVolume();

    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    m_Hit->set_px(1, postPoint->GetMomentum().x() / GeV);
    m_Hit->set_py(1, postPoint->GetMomentum().y() / GeV);
    m_Hit->set_pz(1, postPoint->GetMomentum().z() / GeV);

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
    //sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (!hasMotherSubsystem() && (m_Hit->get_z(1) * cm > m_Zmax || m_Hit->get_z(1) * cm < m_Zmin))
    {
      std::cout << m_Detector->SuperDetector() << std::setprecision(9)
                << " PHG4CylinderSteppingAction: Exit hit z " << m_Hit->get_z(1) * cm
                << " outside acceptance zmin " << m_Zmin
                << ", zmax " << m_Zmax << ", layer: " << layer_id << std::endl;
    }
    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
    }
    else
    {
      if (!m_BlackHoleFlag)
      {
        double eion = edep - aStep->GetNonIonizingEnergyDeposit() / GeV;
        m_Hit->set_eion(m_Hit->get_eion() + eion);
      }
    }
    if (m_SaveLightYieldFlag)
    {
      double light_yield = GetVisibleEnergyDeposition(aStep);
      m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
    }
    if (edep > 0 || m_SaveAllHitsFlag)
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
    if (postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill ||
        m_UseG4StepsFlag > 0)
    {
      // save only hits with energy deposit (or -1 for geantino) or if save all hits flag is set
      if (m_Hit->get_edep() || m_SaveAllHitsFlag)
      {
        m_HitContainer->AddHit(layer_id, m_Hit);
        if (m_SaveShower)
        {
          m_SaveShower->add_g4hit_id(m_HitContainer->GetID(), m_Hit->get_hit_id());
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
void PHG4CylinderSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  // Node Name is passed down from PHG4CylinderSubsystem
  //now look for the map and grab a pointer to it.
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);

  // if we do not find the node we need to scream.
  if (!m_HitContainer && !m_BlackHoleFlag)
  {
    std::cout << "PHG4CylinderSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
  }
}

bool PHG4CylinderSteppingAction::hasMotherSubsystem() const
{
  if (m_Subsystem->GetMotherSubsystem())
  {
    return true;
  }
  return false;
}
