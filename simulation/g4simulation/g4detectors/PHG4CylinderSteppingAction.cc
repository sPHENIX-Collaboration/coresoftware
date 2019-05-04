#include "PHG4CylinderSteppingAction.h"
#include "PHG4CylinderDetector.h"
#include "PHG4StepStatusDecode.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4UserLimits.hh>

#include <boost/io/ios_state.hpp>

#include <iomanip>
#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4CylinderSteppingAction::PHG4CylinderSteppingAction(PHG4CylinderDetector* detector, const PHParameters* parameters)
  : m_Detector(detector)
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
  SetName(m_Detector->GetName());
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
    if ((!isfinite(m_Tmin) && !isfinite(m_Tmax)) ||
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
        aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
    {
      geantino = true;
    }
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();
    //        cout << "time prepoint: " << prePoint->GetGlobalTime()/ns << endl;
    //        cout << "time postpoint: " << postPoint->GetGlobalTime()/ns << endl;
    //        cout << "kinetic energy: " <<  aTrack->GetKineticEnergy()/GeV << endl;
    //       G4ParticleDefinition* def = aTrack->GetDefinition();
    //       cout << "Particle: " << def->GetParticleName() << endl;
    int prepointstatus = prePoint->GetStepStatus();
    if (prepointstatus == fGeomBoundary ||
        prepointstatus == fUndefined ||
        (prepointstatus == fPostStepDoItProc && m_SavePostStepStatus == fGeomBoundary) ||
        m_UseG4StepsFlag > 0)
    {
      if (prepointstatus == fPostStepDoItProc && m_SavePostStepStatus == fGeomBoundary)
      {
        cout << GetName() << ": New Hit for  " << endl;
        cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
             << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
             << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
             << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << endl;
        cout << "last track: " << m_SaveTrackId
             << ", current trackid: " << aTrack->GetTrackID() << endl;
        cout << "phys pre vol: " << volume->GetName()
             << " post vol : " << touchpost->GetVolume()->GetName() << endl;
        cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
             << " previous phys post vol: " << m_SaveVolPost->GetName() << endl;
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

      if (m_Hit->get_z(0) * cm > m_Zmax || m_Hit->get_z(0) * cm < m_Zmin)
      {
        boost::io::ios_precision_saver ips(cout);
        cout << m_Detector->SuperDetector() << std::setprecision(9)
             << "PHG4CylinderSteppingAction: Entry hit z " << m_Hit->get_z(0) * cm
             << " outside acceptance,  zmin " << m_Zmin
             << ", zmax " << m_Zmax << ", layer: " << layer_id << endl;
      }
    }
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
    // some sanity checks for inconsistencies
    // check if this hit was created, if not print out last post step status
    if (!m_Hit || !isfinite(m_Hit->get_x(0)))
    {
      cout << GetName() << ": hit was not created" << endl;
      cout << "prestep status: " << PHG4StepStatusDecode::GetStepStatus(prePoint->GetStepStatus())
           << ", poststep status: " << PHG4StepStatusDecode::GetStepStatus(postPoint->GetStepStatus())
           << ", last pre step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePreStepStatus)
           << ", last post step status: " << PHG4StepStatusDecode::GetStepStatus(m_SavePostStepStatus) << endl;
      cout << "last track: " << m_SaveTrackId
           << ", current trackid: " << aTrack->GetTrackID() << endl;
      cout << "phys pre vol: " << volume->GetName()
           << " post vol : " << touchpost->GetVolume()->GetName() << endl;
      cout << " previous phys pre vol: " << m_SaveVolPre->GetName()
           << " previous phys post vol: " << m_SaveVolPost->GetName() << endl;
      exit(1);
    }
    m_SavePostStepStatus = postPoint->GetStepStatus();
    // check if track id matches the initial one when the hit was created
    if (aTrack->GetTrackID() != m_SaveTrackId)
    {
      cout << "hits do not belong to the same track" << endl;
      cout << "saved track: " << m_SaveTrackId
           << ", current trackid: " << aTrack->GetTrackID()
           << endl;
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
    if (m_Hit->get_z(1) * cm > m_Zmax || m_Hit->get_z(1) * cm < m_Zmin)
    {
      cout << m_Detector->SuperDetector() << std::setprecision(9)
           << " PHG4CylinderSteppingAction: Exit hit z " << m_Hit->get_z(1) * cm
           << " outside acceptance zmin " << m_Zmin
           << ", zmax " << m_Zmax << ", layer: " << layer_id << endl;
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
    if (postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill ||
        m_UseG4StepsFlag > 0)
    {
      // save only hits with energy deposit (or -1 for geantino)
      if (m_Hit->get_edep())
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
  string hitnodename;
  if (m_Detector->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + m_Detector->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + m_Detector->GetName();
  }

  //now look for the map and grab a pointer to it.
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());

  // if we do not find the node we need to make it.
  if (!m_HitContainer && !m_BlackHoleFlag)
  {
    std::cout << "PHG4CylinderSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
}
