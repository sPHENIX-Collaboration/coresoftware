#include "PHG4InnerHcalSteppingAction.h"

#include "PHG4InnerHcalDetector.h"
#include "PHG4StepStatusDecode.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <ffamodules/XploadInterface.h>

#include <phool/getClass.h>

#include <TSystem.h>

// Root headers
#include <TFile.h>
#include <TH2.h>
#include <TSystem.h>

#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>      // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>  // for G4TouchableHandle
#include <Geant4/G4Track.hh>            // for G4Track
#include <Geant4/G4TrackStatus.hh>      // for fStopAndKill
#include <Geant4/G4TransportationManager.hh>
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <cmath>  // for isfinite
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>   // for operator<<, operator+
#include <utility>  // for pair

class PHCompositeNode;

//____________________________________________________________________________..
PHG4InnerHcalSteppingAction::PHG4InnerHcalSteppingAction(PHG4InnerHcalDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_IsActive(m_Params->get_int_param("active"))
  , m_IsBlackHole(m_Params->get_int_param("blackhole"))
  , m_LightScintModelFlag(m_Params->get_int_param("light_scint_model"))
{
  SetLightCorrection(m_Params->get_double_param("light_balance_inner_radius") * cm,
                     m_Params->get_double_param("light_balance_inner_corr"),
                     m_Params->get_double_param("light_balance_outer_radius") * cm,
                     m_Params->get_double_param("light_balance_outer_corr"));
}

PHG4InnerHcalSteppingAction::~PHG4InnerHcalSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
  // since we have a copy in memory of this one - we need to delete it
  delete m_MapCorrHist;
}

//____________________________________________________________________________..
int PHG4InnerHcalSteppingAction::Init()
{
  if (m_LightScintModelFlag)
  {
    std::string ihcalmapname(m_Params->get_string_param("MapFileName"));
    if (ihcalmapname.empty())
    {
      return 0;
    }
    std::string url = XploadInterface::instance()->getUrl("OLD_INNER_HCAL_TILEMAP", ihcalmapname);
    if (!std::filesystem::exists(url))
    {
      std::cout << PHWHERE << " Could not locate " << url << std::endl;
      std::cout << "use empty filename to ignore mapfile" << std::endl;
      gSystem->Exit(1);
    }
    TFile* file = TFile::Open(url.c_str());
    file->GetObject(m_Params->get_string_param("MapHistoName").c_str(), m_MapCorrHist);
    if (!m_MapCorrHist)
    {
      std::cout << "ERROR: could not find Histogram " << m_Params->get_string_param("MapHistoName") << " in " << m_Params->get_string_param("MapFileName") << std::endl;
      gSystem->Exit(1);
      exit(1);  // make code checkers which do not know gSystem->Exit() happy
    }
    m_MapCorrHist->SetDirectory(nullptr);  // rootism: this needs to be set otherwise histo vanished when closing the file
    file->Close();
    delete file;
  }
  return 0;
}

//____________________________________________________________________________..
bool PHG4InnerHcalSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  // m_Detector->IsInInnerHcal(volume)
  // returns
  //  0 is outside of InnerHcal
  //  1 is inside scintillator
  // -1 is steel absorber

  int whichactive = m_Detector->IsInInnerHcal(volume);

  if (!whichactive)
  {
    return false;
  }
  int layer_id = -1;
  int tower_id = -1;
  if (whichactive > 0)  // scintillator
  {
    std::pair<int, int> layer_tower = m_Detector->GetLayerTowerId(volume);
    layer_id = layer_tower.first;
    tower_id = layer_tower.second;
    // std::cout << "name " << volume->GetName() << ", mid: " << layer_id
    //  	   << ", twr: " << tower_id << std::endl;
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
  if (m_IsBlackHole)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  // make sure we are in a volume
  if (m_IsActive)
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
    //       std::cout << "track id " << aTrack->GetTrackID() << std::endl;
    //       std::cout << "time prepoint: " << prePoint->GetGlobalTime() << std::endl;
    //       std::cout << "time postpoint: " << postPoint->GetGlobalTime() << std::endl;
    switch (prePoint->GetStepStatus())
    {
    case fPostStepDoItProc:
      if (m_SavePostStepStatus != fGeomBoundary)
      {
        break;
      }
      else
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
      [[fallthrough]];
    case fGeomBoundary:
    case fUndefined:
      // if previous hit was saved, hit pointer was set to nullptr
      // and we have to make a new one
      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }
      //here we set the entrance values in cm
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);
      // time in ns
      m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
      //set and save the track ID
      m_Hit->set_trkid(aTrack->GetTrackID());
      m_SaveTrackId = aTrack->GetTrackID();
      //set the initial energy deposit
      m_Hit->set_edep(0);
      if (whichactive > 0)  // return of IsInInnerHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
      {
        m_Hit->set_scint_id(tower_id);  // the slat id
        m_Hit->set_eion(0);             // only implemented for v5 otherwise empty
        m_Hit->set_raw_light_yield(0);  //  for scintillator only, initialize light yields
        m_Hit->set_light_yield(0);      // for scintillator only, initialize light yields
        // Now save the container we want to add this hit to
        m_SaveHitContainer = m_Hits;
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
    // some sanity checks for inconsistencies
    // check if this hit was created, if not print out last post step status
    if (!m_Hit || !isfinite(m_Hit->get_x(0)))
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
      gSystem->Exit(1);
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

      gSystem->Exit(1);
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

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);

    //sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (whichactive > 0)  // return of IsInInnerHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
    {
      m_Hit->set_eion(m_Hit->get_eion() + eion);
      light_yield = eion;
      if (m_LightScintModelFlag)
      {
        light_yield = GetVisibleEnergyDeposition(aStep);                         // for scintillator only, calculate light yields
        m_Hit->set_raw_light_yield(m_Hit->get_raw_light_yield() + light_yield);  // save raw Birks light yield
        if (m_MapCorrHist)
        {
          const G4TouchableHandle& theTouchable = prePoint->GetTouchableHandle();
          const G4ThreeVector& worldPosition = postPoint->GetPosition();
          G4ThreeVector localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
          float lx = (localPosition.x() / cm);
          float lz = fabs(localPosition.z() / cm);
          //adjust to tilemap coordinates
          int lcz = (int) (5.0 * lz) + 1;
          int lcx = (int) (5.0 * (lx + 12.1)) + 1;

          if ((lcx >= 1) && (lcx <= m_MapCorrHist->GetNbinsY()) &&
              (lcz >= 1) && (lcz <= m_MapCorrHist->GetNbinsX()))
          {
            light_yield *= (double) (m_MapCorrHist->GetBinContent(lcz, lcx));
          }
          else
          {
            light_yield = 0.0;
          }
        }
        else
        {
          // old correction (linear ligh yield dependence along r), never tested
          light_yield = light_yield * GetLightCorrection(postPoint->GetPosition().x(), postPoint->GetPosition().y());
        }
      }
      m_Hit->set_light_yield(m_Hit->get_light_yield() + light_yield);
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
      if (m_Hit->get_edep() != 0)
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
}

//____________________________________________________________________________..
void PHG4InnerHcalSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  std::string hitnodename;
  std::string absorbernodename;
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
  m_Hits = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  m_AbsorberHits = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename);

  // if we do not find the node it's messed up.
  if (!m_Hits)
  {
    std::cout << "PHG4InnerHcalSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!m_AbsorberHits)
  {
    if (Verbosity() > 1)
    {
      std::cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << absorbernodename << std::endl;
    }
  }
}
