#include "PHG4Prototype2OuterHcalSteppingAction.h"
#include "PHG4Prototype2OuterHcalDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4Prototype2OuterHcalSteppingAction::PHG4Prototype2OuterHcalSteppingAction(PHG4Prototype2OuterHcalDetector* detector, const PHParameters* parameters)
  : m_Detector(detector)
  , m_Hits(nullptr)
  , m_AbsorberHits(nullptr)
  , m_Hit(nullptr)
  , m_Params(parameters)
  , m_SaveHitContainer(nullptr)
  , m_SaveShower(nullptr)
  , m_AbsorberTruthFlag(m_Params->get_int_param("absorbertruth"))
  , m_IsActiveFlag(m_Params->get_int_param("active"))
  , m_IsBlackHoleFlag(m_Params->get_int_param("blackhole"))
  , m_LightScintModelFlag(m_Params->get_int_param("light_scint_model"))
  , m_LightBalanceInnerCorr(m_Params->get_double_param("light_balance_inner_corr"))
  , m_LightBalanceInnerRadius(m_Params->get_double_param("light_balance_inner_radius") * cm)
  , m_LightBalanceOuterCorr(m_Params->get_double_param("light_balance_outer_corr"))
  , m_LightBalanceOuterRadius(m_Params->get_double_param("light_balance_outer_radius") * cm)
{
}

PHG4Prototype2OuterHcalSteppingAction::~PHG4Prototype2OuterHcalSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4Prototype2OuterHcalSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  // m_Detector->IsInPrototype2OuterHcal(volume)
  // returns
  //  0 is outside of Prototype2OuterHcal
  //  1 is inside scintillator
  // -1 is steel absorber

  int whichactive = m_Detector->IsInPrototype2OuterHcal(volume);

  if (!whichactive)
  {
    return false;
  }
  int row_id = -1;
  int slat_id = -1;
  if (whichactive > 0)  // scintillator
  {
    G4VPhysicalVolume* mothervolume = touch->GetVolume(1);
    slat_id = volume->GetCopyNo();
    row_id = m_Detector->get_scinti_row_id(mothervolume->GetName());
     // cout << "mother volume: " <<  mothervolume->GetName()
     //      << ", volume name " << volume->GetName() << ", row: " << row_id
     //  	   << ", column: " << slat_id 
  }
  else
  {
    row_id = m_Detector->get_steel_plate_id(volume->GetName());  // steel plate id
  }
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (m_IsBlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }
  int layer_id = m_Detector->get_Layer();

  // make sure we are in a volume
  if (m_IsActiveFlag)
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
    //       cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
    //       cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:
      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }
      m_Hit->set_row(row_id);  // this is the row
      if (whichactive > 0)   // only for scintillators
      {
        m_Hit->set_scint_id(slat_id);  // the slat id in the mother volume (or steel plate id), the column
      }
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
        }
      }

      //set the initial energy deposit
      m_Hit->set_edep(0);

      m_Hit->set_hit_type(0);
      if ((aTrack->GetParticleDefinition()->GetParticleName().find("e+") != string::npos) ||
          (aTrack->GetParticleDefinition()->GetParticleName().find("e-") != string::npos))
        m_Hit->set_hit_type(1);

      // here we do things which are different between scintillator and absorber hits
      if (whichactive > 0)  // return of IsInPrototype2OuterHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
      {
        m_SaveHitContainer = m_Hits;
        m_Hit->set_light_yield(0);  // for scintillator only, initialize light yields
        m_Hit->set_eion(0);
      }
      else
      {
        m_SaveHitContainer = m_AbsorberHits;
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

    if (whichactive > 0)  // return of IsInPrototype2OuterHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
    {
      m_Hit->set_eion(m_Hit->get_eion() + eion);
      light_yield = eion;
      if (m_LightScintModelFlag)
      {
        light_yield = GetVisibleEnergyDeposition(aStep);  // for scintillator only, calculate light yields
      }
      if (isfinite(m_LightBalanceOuterRadius) &&
          isfinite(m_LightBalanceInnerRadius) &&
          isfinite(m_LightBalanceOuterCorr) &&
          isfinite(m_LightBalanceInnerCorr))
      {
        double r = sqrt(postPoint->GetPosition().x() * postPoint->GetPosition().x() + postPoint->GetPosition().y() * postPoint->GetPosition().y());
        double cor = GetLightCorrection(r);
        light_yield = light_yield * cor;
      }
      m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
    }

    //sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      if (whichactive > 0)
      {
        m_Hit->set_light_yield(-1);
        m_Hit->set_eion(-1);
      }
    }
    if (edep > 0 && (whichactive > 0 || m_AbsorberTruthFlag > 0))
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
    //       m_Hit->identify();
    // return true to indicate the hit was used
    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4Prototype2OuterHcalSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  m_Hits = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  m_AbsorberHits = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!m_Hits)
  {
    std::cout << "PHG4Prototype2OuterHcalSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!m_AbsorberHits)
  {
    if (Verbosity() > 1)
    {
      cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}

double
PHG4Prototype2OuterHcalSteppingAction::GetLightCorrection(const double r) const
{
  double m = (m_LightBalanceOuterCorr - m_LightBalanceInnerCorr) / (m_LightBalanceOuterRadius - m_LightBalanceInnerRadius);
  double b = m_LightBalanceInnerCorr - m * m_LightBalanceInnerRadius;
  double value = m * r + b;
  if (value > 1.0) return 1.0;
  if (value < 0.0) return 0.0;

  return value;
}
