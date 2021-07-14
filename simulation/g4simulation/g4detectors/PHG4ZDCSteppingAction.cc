#include "PHG4ZDCSteppingAction.h"
#include "PHG4ZDCDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4IonisParamMat.hh>  // for G4IonisParamMat
#include <Geant4/G4Material.hh>       // for G4Material
#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <TSystem.h>

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

#include <iostream>
#include <string>  // for basic_string, operator+

class PHCompositeNode;

using namespace std;

//____________________________________________________________________________..
PHG4ZDCSteppingAction::PHG4ZDCSteppingAction(PHG4ZDCDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_SignalHitContainer(nullptr)
  , m_AbsorberHitContainer(nullptr)
  , m_Params(parameters)
  , m_CurrentHitContainer(nullptr)
  , m_Hit(nullptr)
  , m_CurrentShower(nullptr)
  , m_IsActiveFlag(m_Params->get_int_param("active"))
  , absorbertruth(m_Params->get_int_param("absorberactive"))
  , light_scint_model(1)
  , m_IsBlackHole(m_Params->get_int_param("blackhole"))
{
}

PHG4ZDCSteppingAction::~PHG4ZDCSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4ZDCSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch->GetVolume();
 
  // m_Detector->IsInZDC(volume)
  // returns
  //  0 is outside of ZDC
  //  1 is inside scintillator
  // -1 is inside absorber (dead material)

  int whichactive = m_Detector->IsInZDC(volume);
 
  if (!whichactive)
  {
    return false;
  }

  int layer_id = m_Detector->get_Layer();
  int idx_k = -1;
  int idx_j = -1;

  if (whichactive > 0)  // in scintillator
  {
    /* Find indices of scintillator / tower containing this step */
    FindIndex(touch, idx_j, idx_k);
  }
  

  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;

  /* Get pointer to associated Geant4 track */
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (m_IsBlackHole)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  /* Make sure we are in a volume */
  if (m_IsActiveFlag)
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
      
      /* Set hit location (tower index) */
      m_Hit->set_index_k(idx_k);
      m_Hit->set_index_j(idx_j);

      /* Set hit location (space point) */
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

      /* Set hit time */
      m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);

      //set the track ID
      m_Hit->set_trkid(aTrack->GetTrackID());
      /* set intial energy deposit */
      m_Hit->set_edep(0);
      

      /* Now add the hit to the hit collection */
      // here we do things which are different between scintillator and absorber hits
      if (whichactive > 0)
      {
        m_CurrentHitContainer = m_SignalHitContainer;
	m_Hit->set_eion(0);
        m_Hit->set_light_yield(0);  // for scintillator only, initialize light yields
      }
      else
      {
        m_CurrentHitContainer = m_AbsorberHitContainer;
      }
      // here we set what is common for scintillator and absorber hits
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          m_Hit->set_trkid(pp->GetUserTrackId());
          m_Hit->set_shower_id(pp->GetShower()->get_id());
          m_CurrentShower = pp->GetShower();
        }
      }
      break;
    default:
      break;
    }

    if (whichactive > 0)
    {
      light_yield = eion;
      light_yield = GetVisibleEnergyDeposition(aStep);  // for scintillator only, calculate light yields
      static bool once = true;
      if (once && edep > 0)
		 {
          once = false;
	  
          if (Verbosity() > 0)
	    {
	      cout << "PHG4ZDCSteppingAction::UserSteppingAction::"
		//
		   << m_Detector->GetName() << " - "
		   << " use scintillating light model at each Geant4 steps. "
		   << "First step: "
		   << "Material = "
		   << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetName()
		   << ", "
		   << "Birk Constant = "
		   << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetIonisation()->GetBirksConstant()
		   << ","
		   << "edep = " << edep << ", "
		   << "eion = " << eion
		   << ", "
		   << "light_yield = " << light_yield << endl;
	    }
        }
    }
    
    
    /* Update exit values- will be overwritten with every step until
       * we leave the volume or the particle ceases to exist */
    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    /* sum up the energy to get total deposited */
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    m_Hit->set_eion(m_Hit->get_eion() + eion);
    if (whichactive > 0)
    {
      m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
    }

    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      m_Hit->set_eion(-1);
      if (whichactive > 0)
      {
        m_Hit->set_light_yield(-1);
      }
    }
    if (edep > 0 && (whichactive > 0 || absorbertruth > 0))
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
        m_CurrentHitContainer->AddHit(layer_id, m_Hit);
        if (m_CurrentShower)
        {
          m_CurrentShower->add_g4hit_id(m_CurrentHitContainer->GetID(), m_Hit->get_hit_id());
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
void PHG4ZDCSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  m_SignalHitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!m_SignalHitContainer)
  {
    std::cout << "PHG4ZDCSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
    gSystem->Exit(1);
  }
  // this is perfectly fine if absorber hits are disabled
  if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4ZDCSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}

//getting index using copyno
int PHG4ZDCSteppingAction::FindIndex(G4TouchableHandle& touch, int& j, int& k)
{
 

  G4VPhysicalVolume* envelope = touch->GetVolume(2);  //Get the world
  G4VPhysicalVolume* plate = touch->GetVolume(1);  //Get the fiber plate
 
  j = envelope->GetCopyNo();
  k = (plate->GetCopyNo()) / 27;

  return 0;
}

