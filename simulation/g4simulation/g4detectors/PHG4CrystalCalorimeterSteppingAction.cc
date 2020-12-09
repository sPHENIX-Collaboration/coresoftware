#include "PHG4CrystalCalorimeterSteppingAction.h"

#include "PHG4CrystalCalorimeterDefs.h"
#include "PHG4CrystalCalorimeterDetector.h"

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
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation


#include <TSystem.h>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <string>  // for basic_string, string

class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
PHG4CrystalCalorimeterSteppingAction::PHG4CrystalCalorimeterSteppingAction(PHG4CrystalCalorimeterDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_ActiveFlag(parameters->get_int_param("active"))
  , m_BlackHoleFlag(parameters->get_int_param("blackhole"))
{
}

PHG4CrystalCalorimeterSteppingAction::~PHG4CrystalCalorimeterSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4CrystalCalorimeterSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch->GetVolume();

  // m_Detector->IsInCrystalCalorimeter(volume)
  // returns
  //  0 is outside of Crystal Calorimeter
  //  1 is inside scintillator scrystal
  // -1 is absorber (dead material)

  int whichactive = m_Detector->IsInCrystalCalorimeter(volume);

  if (!whichactive)
  {
    return false;
  }

  int layer_id = m_Detector->get_DetectorId();
  int idx_j = -1;
  int idx_k = -1;

  if (whichactive > 0)  // in crystal
  {
    // Find indices of crystal containing this step. The indices are
    // coded into the copy number of the physical volume and we extract them from there
    if (whichactive == PHG4CrystalCalorimeterDefs::CaloType::projective)
    {
      int j[3];
      int k[3];
      for (int i=0;i<3; i++)
      {
	unsigned int icopy = touch->GetVolume(i)->GetCopyNo();
	j[i] =  icopy >> 16;
	k[i] = icopy & 0xFFFF;
      }
      idx_j = j[0] + j[1]*2 + j[2]*4;
      idx_k = k[0] + k[1]*2 + k[2]*4;
    }
    if (whichactive == PHG4CrystalCalorimeterDefs::CaloType::nonprojective)
    {
      unsigned int icopy = touch->GetVolume(1)->GetCopyNo();
      idx_j = icopy >> 16;
      idx_k = icopy & 0xFFFF;
    }
  }

  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
//  G4double light_yield = 0;

  /* Get pointer to associated Geant4 track */
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (m_BlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  /* Make sure we are in a volume */
  if (m_ActiveFlag)
  {
    /* Check if particle is 'geantino' */
    bool geantino = false;
    if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
        aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
    {
      geantino = true;
    }

    /* Get Geant4 pre- and post-step points */
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();

    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:
      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }


      /* Set hit location (space point)*/
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

      /* Set hit time */
      m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);

      /* Set the track ID */
      m_Hit->set_trkid(aTrack->GetTrackID());

      /* Set intial energy deposit */
      m_Hit->set_edep(0);

      /* Now add the hit to the hit collection */
      // here we do things which are different between scintillator and absorber hits
      if (whichactive > 0)
      {
        m_SaveHitContainer = m_HitContainer;
        m_Hit->set_light_yield(0);  // for scintillator only, initialize light yields
        m_Hit->set_eion(0);
        /* Set hit location (tower index) */
        m_Hit->set_index_j(idx_j);
        m_Hit->set_index_k(idx_k);
      }
      else
      {
        m_SaveHitContainer = m_AbsorberHitContainer;
      }

      // here we set what is common for scintillator and absorber hits
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          m_Hit->set_trkid(pp->GetUserTrackId());
          m_Hit->set_shower_id(pp->GetShower()->get_id());
          m_SaveShower = pp->GetShower();
        }
      }
      break;
    default:
      break;
    }

    /* Update exit values- will be overwritten with every step until
       * we leave the volume or the particle ceases to exist */
    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    /* sum up the energy to get total deposited */
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (whichactive > 0)
    {
      m_Hit->set_eion(m_Hit->get_eion() + eion);
      m_Hit->set_light_yield(m_Hit->get_light_yield() + eion);
    }

    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      if (whichactive > 0)
      {
        m_Hit->set_eion(-1);
        m_Hit->set_light_yield(-1);
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
    // (not sure if this will ever be the case)
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
    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4CrystalCalorimeterSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename);

  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4CrystalCalorimeterSteppingAction::SetInterfacePointers - unable to find " << hitnodename << std::endl;
    gSystem->Exit(1);
  }
  if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4CrystalCalorimeterSteppingAction::SetInterfacePointers - unable to find " << absorbernodename << endl;
    }
  }
}
