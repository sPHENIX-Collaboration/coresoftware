// local headers in quotes (that is important when using include subdirs!)
#include "PHG4OuterHcalSteppingAction.h"
#include "PHG4HcalDefs.h"
#include "PHG4OuterHcalDetector.h"
#include "PHG4StepStatusDecode.h"

// our own headers in alphabetical order

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllServer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

// Geant4 headers

#include <Geant4/G4Field.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4PropagatorInField.hh>
#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4TransportationManager.hh>

// Root headers
#include <TH2F.h>
#include <TSystem.h>

// boost headers
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

// finally system headers
#include <cassert>
#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4OuterHcalSteppingAction::PHG4OuterHcalSteppingAction(PHG4OuterHcalDetector* detector, const PHParameters* parameters)
  : detector_(detector)
  , hits_(nullptr)
  , absorberhits_(nullptr)
  , hit(nullptr)
  , params(parameters)
  , savehitcontainer(nullptr)
  , saveshower(nullptr)
  , savevolpre(nullptr)
  , savevolpost(nullptr)
  , savetrackid(-1)
  , saveprestepstatus(-1)
  , savepoststepstatus(-1)
  , enable_field_checker(0)
  , absorbertruth(params->get_int_param("absorbertruth"))
  , IsActive(params->get_int_param("active"))
  , IsBlackHole(params->get_int_param("blackhole"))
  , n_scinti_plates(params->get_int_param(PHG4HcalDefs::scipertwr) * params->get_int_param("n_towers"))
  , light_scint_model(params->get_int_param("light_scint_model"))
  , light_balance_inner_corr(params->get_double_param("light_balance_inner_corr"))
  , light_balance_inner_radius(params->get_double_param("light_balance_inner_radius") * cm)
  , light_balance_outer_corr(params->get_double_param("light_balance_outer_corr"))
  , light_balance_outer_radius(params->get_double_param("light_balance_outer_radius") * cm)
{
  name = detector_->GetName();
}

PHG4OuterHcalSteppingAction::~PHG4OuterHcalSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

int PHG4OuterHcalSteppingAction::Init()
{
  enable_field_checker = GetIntOpt("FieldChecker");
  return 0;
}

//____________________________________________________________________________..
bool PHG4OuterHcalSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  // detector_->IsInOuterHcal(volume)
  // returns
  //  0 is outside of OuterHcal
  //  1 is inside scintillator
  // -1 is steel absorber (if absorber set to active)

  int whichactive = detector_->IsInOuterHcal(volume);

  if (!whichactive)
  {
    return false;
  }

  if (enable_field_checker)
  {
    FieldChecker(aStep);
  }

  int layer_id = -1;
  int tower_id = -1;
  if (whichactive > 0)  // scintillator
  {
    // G4AssemblyVolumes naming convention:
    //     av_WWW_impr_XXX_YYY_ZZZ
    // where:

    //     WWW - assembly volume instance number
    //     XXX - assembly volume imprint number
    //     YYY - the name of the placed logical volume
    //     ZZZ - the logical volume index inside the assembly volume
    // e.g. av_1_impr_82_HcalOuterScinti_11_pv_11
    // 82 the number of the scintillator mother volume
    // HcalOuterScinti_11: name of scintillator slat
    // 11: number of scintillator slat logical volume
    // use boost tokenizer to separate the _, then take value
    // after "impr" for mother volume and after "pv" for scintillator slat
    // use boost lexical cast for string -> int conversion
    boost::char_separator<char> sep("_");
    boost::tokenizer<boost::char_separator<char> > tok(volume->GetName(), sep);
    boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
    for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
    {
      if (*tokeniter == "impr")
      {
        ++tokeniter;
        if (tokeniter != tok.end())
        {
          layer_id = boost::lexical_cast<int>(*tokeniter);
          // check detector description, for assemblyvolumes it is not possible
          // to give the first volume id=0, so they go from id=1 to id=n.
          // I am not going to start with fortran again - our indices start
          // at zero, id=0 to id=n-1. So subtract one here
          layer_id--;
          if (layer_id < 0 || layer_id >= n_scinti_plates)
          {
            cout << "invalid scintillator row " << layer_id
                 << ", valid range 0 < row < " << n_scinti_plates << endl;
            gSystem->Exit(1);
          }
        }
        else
        {
          cout << PHWHERE << " Error parsing " << volume->GetName()
               << " for mother volume number " << endl;
          gSystem->Exit(1);
        }
      }
      else if (*tokeniter == "pv")
      {
        ++tokeniter;
        if (tokeniter != tok.end())
        {
          tower_id = boost::lexical_cast<int>(*tokeniter);
        }
        else
        {
          cout << PHWHERE << " Error parsing " << volume->GetName()
               << " for mother scinti slat id " << endl;
          gSystem->Exit(1);
        }
      }
    }
  }
  else
  {
    layer_id = touch->GetCopyNumber();  // steel plate id
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
    case fPostStepDoItProc:
      if (savepoststepstatus != fGeomBoundary)
      {
        break;
      }
      else
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
    case fGeomBoundary:
    case fUndefined:
      if (!hit)
      {
        hit = new PHG4Hitv1();
      }
      //here we set the entrance values in cm
      hit->set_x(0, prePoint->GetPosition().x() / cm);
      hit->set_y(0, prePoint->GetPosition().y() / cm);
      hit->set_z(0, prePoint->GetPosition().z() / cm);
      // time in ns
      hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
      //set the track ID
      hit->set_trkid(aTrack->GetTrackID());
      savetrackid = aTrack->GetTrackID();
      //set the initial energy deposit
      hit->set_edep(0);
      if (whichactive > 0)  // return of IsInOuterHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
      {
        hit->set_scint_id(tower_id);  // the slat id
        hit->set_eion(0);
        hit->set_light_yield(0);  //  for scintillator only, initialize light yields
        // Now save the container we want to add this hit to
        savehitcontainer = hits_;
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
    // some sanity checks for inconsistencies
    // check if this hit was created, if not print out last post step status
    if (!hit || !isfinite(hit->get_x(0)))
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
    savepoststepstatus = postPoint->GetStepStatus();
    // check if track id matches the initial one when the hit was created
    if (aTrack->GetTrackID() != savetrackid)
    {
      cout << GetName() << ": hits do not belong to the same track" << endl;
      cout << "saved track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID()
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
    hit->set_x(1, postPoint->GetPosition().x() / cm);
    hit->set_y(1, postPoint->GetPosition().y() / cm);
    hit->set_z(1, postPoint->GetPosition().z() / cm);

    hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    if (whichactive > 0)
    {
      if (light_scint_model)
      {
        light_yield = GetVisibleEnergyDeposition(aStep);
      }
      else
      {
        light_yield = eion;
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
    }

    //sum up the energy to get total deposited
    hit->set_edep(hit->get_edep() + edep);
    if (whichactive > 0)
    {
      hit->set_eion(hit->get_eion() + eion);
      hit->set_light_yield(hit->get_light_yield() + light_yield);
    }
    if (geantino)
    {
      hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      hit->set_eion(-1);
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
    // return true to indicate the hit was used
    return true;
  }
  else
  {
    return false;
  }
}

//____________________________________________________________________________..
void PHG4OuterHcalSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  absorberhits_ = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!hits_)
  {
    std::cout << "PHG4OuterHcalSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
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
PHG4OuterHcalSteppingAction::GetLightCorrection(const double r) const
{
  double m = (light_balance_outer_corr - light_balance_inner_corr) / (light_balance_outer_radius - light_balance_inner_radius);
  double b = light_balance_inner_corr - m * light_balance_inner_radius;
  double value = m * r + b;
  if (value > 1.0) return 1.0;
  if (value < 0.0) return 0.0;

  return value;
}

void PHG4OuterHcalSteppingAction::FieldChecker(const G4Step* aStep)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  assert(se);

  static const string h_field_name = "hOuterHcalField";

  if (not se->isHistoRegistered(h_field_name))
  {
    TH2F* h = new TH2F(h_field_name.c_str(), "Magnetic field (Tesla) in HCal;X (cm);Y (cm)", 2400,
                       -300, 300, 2400, -300, 300);

    se->registerHisto(h, 1);

    cout << "PHG4OuterHcalSteppingAction::FieldChecker - make a histograme to check outer Hcal field map."
         << " Saved to Fun4AllServer Histo with name " << h_field_name << endl;
  }

  TH2F* h = dynamic_cast<TH2F*>(se->getHisto(h_field_name));
  assert(h);

  assert(aStep);
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  assert(touch);
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();
  assert(volume);

  G4ThreeVector globPosition = aStep->GetPreStepPoint()->GetPosition();

  G4double globPosVec[4] = {0};
  G4double FieldValueVec[6] = {0};

  globPosVec[0] = globPosition.x();
  globPosVec[1] = globPosition.y();
  globPosVec[2] = globPosition.z();
  globPosVec[3] = aStep->GetPreStepPoint()->GetGlobalTime();

  const Int_t binx = h->GetXaxis()->FindBin(globPosVec[0] / cm);
  const Int_t biny = h->GetYaxis()->FindBin(globPosVec[1] / cm);

  if (h->GetBinContent(binx, binx) == 0)
  {  // only fille unfilled bins

    G4TransportationManager* transportMgr =
        G4TransportationManager::GetTransportationManager();
    assert(transportMgr);

    G4PropagatorInField* fFieldPropagator =
        transportMgr->GetPropagatorInField();
    assert(fFieldPropagator);

    G4FieldManager* fieldMgr = fFieldPropagator->FindAndSetFieldManager(volume);
    assert(fieldMgr);

    const G4Field* pField = fieldMgr->GetDetectorField();
    assert(pField);

    pField->GetFieldValue(globPosVec, FieldValueVec);

    G4ThreeVector FieldValue = G4ThreeVector(FieldValueVec[0],
                                             FieldValueVec[1], FieldValueVec[2]);

    const double B = FieldValue.mag() / tesla;

    h->SetBinContent(binx, biny, B);

    cout << "PHG4OuterHcalSteppingAction::FieldChecker - "
         << "bin " << binx
         << ", " << biny << " := " << B << " Tesla @ x,y = " << globPosVec[0] / cm
         << "," << globPosVec[1] / cm << " cm" << endl;
  }
}
