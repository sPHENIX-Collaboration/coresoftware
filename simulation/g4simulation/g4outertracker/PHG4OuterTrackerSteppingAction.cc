#include "PHG4OuterTrackerSteppingAction.h"

#include "PHG4OuterTrackerDetector.h"

#include <g4detectors/PHG4StepStatusDecode.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>         // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>
#include <phool/phool.h>                       // for PHWHERE

#include <Geant4/G4NavigationHistory.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>               // for G4StepPoint
#include <Geant4/G4StepStatus.hh>              // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>                  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>             // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>         // for G4TouchableHandle
#include <Geant4/G4Track.hh>                   // for G4Track
#include <Geant4/G4TrackStatus.hh>             // for fStopAndKill
#include <Geant4/G4Types.hh>                   // for G4double
#include <Geant4/G4VPhysicalVolume.hh>         // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>              // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>   // for G4VUserTrackInformation

#include <boost/tokenizer.hpp>
// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp>  // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 4 && BOOST_VERSION == 105700)
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma message "ignoring bogus gcc warning in boost header lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#pragma GCC diagnostic warning "-Wuninitialized"
#else
#include <boost/lexical_cast.hpp>
#endif

#include <cmath>                              // for NAN
#include <cstdlib>                            // for exit
#include <iostream>
#include <string>                              // for operator<<, basic_string

class PHCompositeNode;

using namespace std;
//____________________________________________________________________________..
PHG4OuterTrackerSteppingAction::PHG4OuterTrackerSteppingAction(PHG4OuterTrackerDetector* detector,  const int layer)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_HitContainer(nullptr)
  , m_AbsorberhitContainer(nullptr)
  , m_Hit(nullptr)
  , m_SaveShower(nullptr)
{
  //params = _paramsContainer->GetParameters(layer);
  layer_id = layer;
}

//____________________________________________________________________________..
PHG4OuterTrackerSteppingAction::~PHG4OuterTrackerSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4OuterTrackerSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{

  //=======================================================================
  // We want the location of the hit
  // Here we will record in the hit object:
  //   The energy deposited
  //   The entry point and exit point in world coordinates
  //   The entry point and exit point in local (sensor) coordinates
  // The pixel number will be derived later from the entry and exit points in the sensor local coordinates
  //=======================================================================

  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  //  if (Verbosity() > 5)
  {
    // make sure we know where we are!
    G4VPhysicalVolume* vtest = touch->GetVolume();
    if(Verbosity() > 100) cout << "Entering PHG4OuterTrackerSteppingAction::UserSteppingAction, checking for volume " << vtest->GetName() << endl;
  }

  int whichactive = m_Detector->IsInOuterTracker(volume);
  if (!whichactive)
  {
    if(Verbosity() > 100) cout << "   skip this volume, not in " << m_Detector->GetName() << endl;
    return false;
  }

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  int IsBlackHole = m_Detector->IsBlackHole(volume);
  if (IsBlackHole)
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
      aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
  {
    geantino = true;
  }
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  int prepointstatus = prePoint->GetStepStatus();

  if(Verbosity() > 5)
    {
      cout << "track id " << aTrack->GetTrackID() << endl;
      cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
      cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
    }

  if ((whichactive > 0) ||
      prepointstatus == fGeomBoundary ||
      prepointstatus == fUndefined ||
      (prepointstatus == fPostStepDoItProc && savepoststepstatus == fGeomBoundary))
  {
    // this is for debugging weird occurances we have occasionally
    if (prepointstatus == fPostStepDoItProc && savepoststepstatus == fGeomBoundary)
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

    StoreLocalCoordinate(m_Hit, aStep, true, false);

    // time in ns
    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
    //set and save the track ID
    m_Hit->set_trkid(aTrack->GetTrackID());
    savetrackid = aTrack->GetTrackID();
    //set the initial energy deposit
    m_Hit->set_edep(0);
    if (whichactive > 0)  
    {
      m_Hit->set_eion(0);
      // Now save the container we want to add this m_Hit to
      m_SaveHitContainer = m_HitContainer;
    }
    else
    {
      m_SaveHitContainer = m_AbsorberhitContainer;
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
    // some sanity checks for inconsistencies
    // check if this hit was created, if not print out last post step status
    if (!m_Hit || !isfinite(m_Hit->get_x(0)))
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
    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    StoreLocalCoordinate(m_Hit, aStep, false, true);
    
    if (whichactive > 0)
      {
	m_Hit->set_px(1, postPoint->GetMomentum().x() / GeV);
	m_Hit->set_py(1, postPoint->GetMomentum().y() / GeV);
	m_Hit->set_pz(1, postPoint->GetMomentum().z() / GeV);
	m_Hit->set_eion(m_Hit->get_eion() + eion);
      }

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    //sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (whichactive > 0)
    {
      m_Hit->set_eion(m_Hit->get_eion() + eion);
    }
    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4m_Hits get dropped, this way geantinos survive the g4hit compression
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
    if ((whichactive > 0) ||
        postPoint->GetStepStatus() == fGeomBoundary ||
        postPoint->GetStepStatus() == fWorldBoundary ||
        postPoint->GetStepStatus() == fAtRestDoItProc ||
        aTrack->GetTrackStatus() == fStopAndKill)
    {
      // save only hits with energy deposit (or -1 for geantino)
      if (m_Hit->get_edep())
      {
        m_SaveHitContainer->AddHit(layer_id, m_Hit);
        if (m_SaveShower)
        {
          m_SaveShower->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());
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

  return false;
}

//____________________________________________________________________________..
void PHG4OuterTrackerSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  string hitnodename;
  string absorbernodename;
  if (m_Detector->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + m_Detector->SuperDetector();
    absorbernodename = "G4HIT_ABSORBER_" + m_Detector->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + m_Detector->GetName();
    absorbernodename = "G4HIT_ABSORBER_" + m_Detector->GetName();
  }

  //now look for the map and grab a pointer to it.
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  m_AbsorberhitContainer = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4OuterTrackerSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!m_AbsorberhitContainer)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4OuterTrackerSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}
