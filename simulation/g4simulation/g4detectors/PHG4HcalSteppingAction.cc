#include "PHG4HcalSteppingAction.h"
#include "PHG4HcalDetector.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4ParticleDefinition.hh>  // for G4ParticleDefinition
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>              // for G4StepPoint
#include <Geant4/G4StepStatus.hh>             // for fGeomBoundary, fAtRestD...
#include <Geant4/G4String.hh>                 // for G4String
#include <Geant4/G4SystemOfUnits.hh>          // for cm, GeV, nanosecond
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>        // for G4TouchableHandle
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <iostream>
#include <string>  // for char_traits, operator<<

class G4VPhysicalVolume;
class PHCompositeNode;

using namespace std;
//____________________________________________________________________________..
PHG4HcalSteppingAction::PHG4HcalSteppingAction(PHG4HcalDetector* detector)
  : PHG4SteppingAction(detector->GetName())
  , detector_(detector)
{
}

//____________________________________________________________________________..
PHG4HcalSteppingAction::~PHG4HcalSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4HcalSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  // get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;

  const G4Track* aTrack = aStep->GetTrack();

  int layer_id = detector_->get_Layer();
  // make sure we are in a volume
  // IsInCylinderActive returns the number of the scintillator
  // slat which has fired
  int isactive = detector_->IsInCylinderActive(volume);
  if (isactive > PHG4HcalDetector::INACTIVE)
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
    //       cout << "track id " << aTrack->GetTrackID() << endl;
    //        cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
    //        cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:

      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }
      m_Hit->set_layer((unsigned int) layer_id);
      m_Hit->set_scint_id(isactive);  // isactive contains the scintillator slat id
      //here we set the entrance values in cm
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

      // time in ns
      m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
      //set the track ID
      m_Hit->set_trkid(aTrack->GetTrackID());
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          m_Hit->set_trkid(pp->GetUserTrackId());
          m_Hit->set_shower_id(pp->GetShower()->get_id());
          m_SaveShower = pp->GetShower();
        }
      }

      //set the initial energy deposit
      m_Hit->set_edep(0);
      m_Hit->set_eion(0);  // only implemented for v5 otherwise empty
      m_Hit->set_light_yield(0);
      //	  m_Hit->print();
      // Now add the hit
      if (isactive >= 0)  // the slat ids start with zero
      {
        m_SaveHitContainer = m_HitContainer;
        //	      unsigned int shift_layer_id = layer_id << (phg4hitdefs::keybits - 3);
      }
      else
      {
        m_SaveHitContainer = m_AbsorberHits;
      }
      //      m_SaveHitContainer->AddHit(layer_id, m_Hit);

      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          pp->GetShower()->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());
        }
      }
      if (m_Hit->get_z(0) > zmax || m_Hit->get_z(0) < zmin)
      {
        cout << "PHG4HcalSteppingAction: hit outside acceptance, layer: " << layer_id << endl;
        m_Hit->identify();
      }
      break;
    default:
      break;
    }
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    if (isactive >= 0)  // the slat ids start with zero
    {
      if (light_scint_model_)
      {
        light_yield = GetVisibleEnergyDeposition(aStep);
      }
      else
      {
        light_yield = eion;
      }

      if (ValidCorrection())
      {
        light_yield = light_yield * GetLightCorrection(postPoint->GetPosition().x() / cm, (postPoint->GetPosition().y() / cm));
      }
    }

    //sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    m_Hit->set_eion(m_Hit->get_eion() + eion);
    m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
    m_Hit->set_path_length(aTrack->GetTrackLength() / cm);

    if (m_Hit->get_z(1) > zmax || m_Hit->get_z(1) < zmin)
    {
      cout << "PHG4HcalSteppingAction: hit outside acceptance zmin " << zmin << ", zmax " << zmax << " at exit" << endl;
      m_Hit->identify();
    }
    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      m_Hit->set_eion(-1);
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
    //       hit->identify();
    // return true to indicate the hit was used
    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4HcalSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  string hitnodename;
  string absorbernodename;
  if (detector_->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + detector_->SuperDetector();
    absorbernodename = "G4HIT_ABSORBER_" + detector_->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + detector_->GetName();
    absorbernodename = "G4HIT_ABSORBER_" + detector_->GetName();
  }

  //now look for the map and grab a pointer to it.
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  m_AbsorberHits = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename);
  // if we do not find the node we need to make it.
  if (!m_HitContainer)
  {
    std::cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!m_AbsorberHits)
  {
    if (Verbosity() > 0)
    {
      std::cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << absorbernodename << std::endl;
    }
  }
}
