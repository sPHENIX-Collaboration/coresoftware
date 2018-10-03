#include "PHG4Prototype2InnerHcalSteppingAction.h"
#include "PHG4Prototype2InnerHcalDetector.h"

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

#include <boost/foreach.hpp>
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

//#define TESTSINGLESLAT
#ifdef TESTSINGLESLAT
static const int nrow = 4;
static const int nslat = 9;
#endif

using namespace std;
//____________________________________________________________________________..
PHG4Prototype2InnerHcalSteppingAction::PHG4Prototype2InnerHcalSteppingAction(PHG4Prototype2InnerHcalDetector* detector, const PHParameters* parameters)
  : m_Detector(detector)
  , hits_(nullptr)
  , absorberhits_(nullptr)
  , hit(nullptr)
  , params(parameters)
  , savehitcontainer(nullptr)
  , saveshower(nullptr)
  , absorbertruth(params->get_int_param("absorbertruth"))
  , IsActive(params->get_int_param("active"))
  , IsBlackHole(params->get_int_param("blackhole"))
  , light_scint_model(params->get_int_param("light_scint_model"))
  , light_balance_inner_corr(params->get_double_param("light_balance_inner_corr"))
  , light_balance_inner_radius(params->get_double_param("light_balance_inner_radius") * cm)
  , light_balance_outer_corr(params->get_double_param("light_balance_outer_corr"))
  , light_balance_outer_radius(params->get_double_param("light_balance_outer_radius") * cm)
{
}

PHG4Prototype2InnerHcalSteppingAction::~PHG4Prototype2InnerHcalSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4Prototype2InnerHcalSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  // m_Detector->IsInPrototype2InnerHcal(volume)
  // returns
  //  0 is outside of Prototype2InnerHcal
  //  1 is inside scintillator
  // -1 is steel absorber

  int whichactive = m_Detector->IsInPrototype2InnerHcal(volume);

  if (!whichactive)
  {
    return false;
  }
  int row_id = -1;
  int slat_id = -1;
  if (whichactive > 0)  // scintillator
  {
    // first extract the scintillator id (0-3) from the volume name (OuterScinti_0,1,2,3)
    boost::char_separator<char> sep("_");
    boost::tokenizer<boost::char_separator<char> > tok(volume->GetName(), sep);
    boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter = tok.begin();
    ++tokeniter;
    slat_id = boost::lexical_cast<int>(*tokeniter);
    // G4AssemblyVolumes naming convention:
    //     av_WWW_impr_XXX_YYY_ZZZ
    // where:

    //     WWW - assembly volume instance number
    //     XXX - assembly volume imprint number
    //     YYY - the name of the placed logical volume
    //     ZZZ - the logical volume index inside the assembly volume
    // e.g. av_1_impr_2_InnerHcalScintiMother_pv_11
    // 2 the number of the scintillator mother volume
    // InnerHcalScintiMother_11: name of scintillator slat
    // 11: number of scintillator slat logical volume
    // use boost tokenizer to separate the _, then take value
    // after "impr" for mother volume and after "pv" for scintillator slat
    // use boost lexical cast for string -> int conversion
    G4VPhysicalVolume* mothervolume = touch->GetVolume(1);
    boost::tokenizer<boost::char_separator<char> > tokm(mothervolume->GetName(), sep);
    for (tokeniter = tokm.begin(); tokeniter != tokm.end(); ++tokeniter)
    {
      if (*tokeniter == "pv")
      {
        ++tokeniter;
        if (tokeniter != tokm.end())
        {
          row_id = boost::lexical_cast<int>(*tokeniter);
          // from the construction via assembly volumes, the mother id starts
          // at 2 and is incremented by 2 for each new row of slats
          // this maps it back to 0-nslats
          row_id -= 2;
          row_id /= 2;
        }
        else
        {
          cout << PHWHERE << " Error parsing " << mothervolume->GetName()
               << " for mother scinti slat id " << endl;
          exit(1);
        }
        break;
      }
    }
    if (slat_id != volume->GetCopyNo() || row_id != m_Detector->get_scinti_row_id(mothervolume->GetName()))
    {
    cout << "mother volume: " <<  mothervolume->GetName()
          << ", volume name " << volume->GetName() << ", row: " << row_id
	 << ", 2: " << m_Detector->get_scinti_row_id(mothervolume->GetName())
      	   << ", column: " << slat_id 
	   << ", 2: " << volume->GetCopyNo() << endl;
    exit(1);
    }
#ifdef TESTSINGLESLAT
    if (row_id != nrow)
    {
      return false;
    }
    if (slat_id != nslat)
    {
      return false;
    }
#endif
  }
  else
  {
    slat_id = touch->GetCopyNumber();  // steel plate id
#ifdef TESTSINGLESLAT
    return false;
#endif
  }
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (IsBlackHole)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }
  int layer_id = m_Detector->get_Layer();
  // make sure we are in a volume
  if (IsActive)
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
      if (!hit)
      {
        hit = new PHG4Hitv1();
      }
      hit->set_row(row_id);
      if (whichactive > 0)  // only for scintillators
      {
        hit->set_scint_id(slat_id);  // the slat id in the mother volume (or steel plate id), the column
      }
      //here we set the entrance values in cm
      hit->set_x(0, prePoint->GetPosition().x() / cm);
      hit->set_y(0, prePoint->GetPosition().y() / cm);
      hit->set_z(0, prePoint->GetPosition().z() / cm);
      // time in ns
      hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
      //set the track ID
      hit->set_trkid(aTrack->GetTrackID());
      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          hit->set_trkid(pp->GetUserTrackId());
          hit->set_shower_id(pp->GetShower()->get_id());
        }
      }

      //set the initial energy deposit
      hit->set_edep(0);

      hit->set_hit_type(0);
      if ((aTrack->GetParticleDefinition()->GetParticleName().find("e+") != string::npos) ||
          (aTrack->GetParticleDefinition()->GetParticleName().find("e-") != string::npos))
        hit->set_hit_type(1);

      if (whichactive > 0)  // return of IsInPrototype2InnerHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
      {
        savehitcontainer = hits_;
        hit->set_light_yield(0);  // for scintillator only, initialize light yields
        hit->set_eion(0);
      }
      else
      {
        savehitcontainer = absorberhits_;
      }

      if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
      {
        if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
        {
          hit->set_trkid(pp->GetUserTrackId());
          hit->set_shower_id(pp->GetShower()->get_id());
          saveshower = pp->GetShower();
        }
      }
      break;
    default:
      break;
    }
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
    hit->set_x(1, postPoint->GetPosition().x() / cm);
    hit->set_y(1, postPoint->GetPosition().y() / cm);
    hit->set_z(1, postPoint->GetPosition().z() / cm);

    hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    if (whichactive > 0)  // return of IsInPrototype2InnerHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
    {
      hit->set_eion(hit->get_eion() + eion);
      light_yield = eion;
      if (light_scint_model)
      {
        light_yield = GetVisibleEnergyDeposition(aStep);  // for scintillator only, calculate light yields
      }
      if (isfinite(light_balance_outer_radius) &&
          isfinite(light_balance_inner_radius) &&
          isfinite(light_balance_outer_corr) &&
          isfinite(light_balance_inner_corr))
      {
        double r = sqrt(postPoint->GetPosition().x() * postPoint->GetPosition().x() + postPoint->GetPosition().y() * postPoint->GetPosition().y());
        double cor = GetLightCorrection(r);
        light_yield = light_yield * cor;
      }
      hit->set_light_yield(hit->get_light_yield() + light_yield);
    }

    //sum up the energy to get total deposited
    hit->set_edep(hit->get_edep() + edep);
    if (geantino)
    {
      hit->set_edep(-1);    // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      if (whichactive > 0)  // add light yield for scintillators
      {
        hit->set_light_yield(-1);
        hit->set_eion(-1);
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
      if (hit->get_edep())
      {
        savehitcontainer->AddHit(layer_id, hit);
        if (saveshower)
        {
          saveshower->add_g4hit_id(savehitcontainer->GetID(), hit->get_hit_id());
        }
        // ownership has been transferred to container, set to null
        // so we will create a new hit for the next track
        hit = nullptr;
      }
      else
      {
        // if this hit has no energy deposit, just reset it for reuse
        // this means we have to delete it in the dtor. If this was
        // the last hit we processed the memory is still allocated
        hit->Reset();
      }
    }

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
void PHG4Prototype2InnerHcalSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  absorberhits_ = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!hits_)
  {
    std::cout << "PHG4Prototype2InnerHcalSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!absorberhits_)
  {
    if (Verbosity() > 1)
    {
      cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}

double
PHG4Prototype2InnerHcalSteppingAction::GetLightCorrection(const double r) const
{
  double m = (light_balance_outer_corr - light_balance_inner_corr) / (light_balance_outer_radius - light_balance_inner_radius);
  double b = light_balance_inner_corr - m * light_balance_inner_radius;
  double value = m * r + b;
  if (value > 1.0) return 1.0;
  if (value < 0.0) return 0.0;

  return value;
}
