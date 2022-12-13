
// local headers in quotes (that is important when using include subdirs!)
#include "PHG4OHCalSteppingAction.h"

#include "PHG4OHCalDetector.h"

// our own headers in alphabetical order

#include <g4detectors/PHG4HcalDefs.h>
#include <g4detectors/PHG4StepStatusDecode.h>

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllServer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

// Root headers
#include <TAxis.h>  // for TAxis
#include <TFile.h>
#include <TH2.h>
#include <TNamed.h>  // for TNamed
#include <TSystem.h>

// Geant4 headers

#include <Geant4/G4AffineTransform.hh>  // for G4AffineTransform
#include <Geant4/G4Field.hh>
#include <Geant4/G4FieldManager.hh>
#include <Geant4/G4NavigationHistory.hh>   // for G4NavigationHistory
#include <Geant4/G4ParticleDefinition.hh>  // for G4ParticleDefinition
#include <Geant4/G4PropagatorInField.hh>
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

// finally system headers
#include <cassert>
#include <cmath>  // for isfinite, sqrt
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>  // for operator<<, string
#include <tuple>   // for get, tuple

class PHCompositeNode;

//____________________________________________________________________________..
PHG4OHCalSteppingAction::PHG4OHCalSteppingAction(PHG4OHCalDetector* detector, const PHParameters* parameters)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_Params(parameters)
  , m_EnableFieldCheckerFlag(m_Params->get_int_param("field_check"))
  , m_IsActiveFlag(m_Params->get_int_param("active"))
  , m_IsBlackHoleFlag(m_Params->get_int_param("blackhole"))
  , m_NScintiPlates(m_Params->get_int_param(PHG4HcalDefs::scipertwr) * m_Params->get_int_param("n_towers"))
  , m_LightScintModelFlag(m_Params->get_int_param("light_scint_model"))
{
  SetLightCorrection(m_Params->get_double_param("light_balance_inner_radius") * cm,
                     m_Params->get_double_param("light_balance_inner_corr"),
                     m_Params->get_double_param("light_balance_outer_radius") * cm,
                     m_Params->get_double_param("light_balance_outer_corr"));
}

PHG4OHCalSteppingAction::~PHG4OHCalSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;

  // since we have a copy in memory of this one - we need to delete it
  for (auto it : m_MapCorrHist)
  {
    delete it;
  }
  for (int j = 0; j < 4; j++)
  {
    delete m_MapCorrHistChim[j];
  }
}

int PHG4OHCalSteppingAction::Init()
{
  m_EnableFieldCheckerFlag = m_Params->get_int_param("field_check");

  if (m_LightScintModelFlag)
  {
    std::string mappingfilename(m_Params->get_string_param("MapFileName"));
    if (mappingfilename.empty())
    {
      return 0;
    }
    if (!std::filesystem::exists(m_Params->get_string_param("MapFileName")))
    {
      std::cout << PHWHERE << " Could not locate " << m_Params->get_string_param("MapFileName") << std::endl;
      std::cout << "use empty filename to ignore mapfile" << std::endl;
      gSystem->Exit(1);
    }

    TFile* file = TFile::Open(mappingfilename.c_str());
    std::string Tilehist = m_Params->get_string_param("MapHistoName");

    for (int i = 0; i < 24; i++)
    {
      std::string str2 = std::to_string(i);
      Tilehist += str2;
      file->GetObject(Tilehist.c_str(), m_MapCorrHist[i]);
      if (i < 4)
      {
        Tilehist += "_chimney";
      }
      file->GetObject(Tilehist.c_str(), m_MapCorrHistChim[i]);
      Tilehist = m_Params->get_string_param("MapHistoName");

      if ((!m_MapCorrHist[i]) || (!m_MapCorrHistChim[i]))
      {
        std::cout << "ERROR: could not find Histogram " << Tilehist << i << " in " << m_Params->get_string_param("MapFileName") << std::endl;
        gSystem->Exit(1);
      }

      m_MapCorrHist[i]->SetDirectory(nullptr);  // rootism: this needs to be set otherwise histo vanished when closing the file
      m_MapCorrHistChim[i]->SetDirectory(nullptr);
    }

    file->Close();
    delete file;
  }
  return 0;
}

//____________________________________________________________________________..
bool PHG4OHCalSteppingAction::UserSteppingAction(const G4Step* aStep, bool /*was_used*/)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  // m_Detector->IsInOHCal(volume)
  // returns
  //  0 is outside of OHCal
  //  1 is inside scintillator
  // -1 is steel absorber (if absorber set to active)

  int whichactive = m_Detector->IsInOHCal(volume);

  if (!whichactive)
  {
    return false;
  }

  if (m_EnableFieldCheckerFlag)
  {
    FieldChecker(aStep);
  }

  int layer_id = -1;
  int tower_id = -1;
  int sector_id = -1;
  if (whichactive > 0)  // scintillator
  {
    std::tuple<int, int, int> layer_tower = m_Detector->GetRowColumnId(volume);
    sector_id = std::get<0>(layer_tower);  // calorimeter sector (0-31)
    layer_id = std::get<1>(layer_tower);   // calorimeter scintillator row (0-319)
    tower_id = std::get<2>(layer_tower);   // calorimeter scintillator tower in row (0-23)
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
  if (m_IsBlackHoleFlag)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  // make sure we are in a volume
  if (m_IsActiveFlag)
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
        if (m_SavePostStepStatus != fAtRestDoItProc)
        {
          break;
        }
        else
        {
          if (aTrack->GetTrackID() == m_SaveTrackId)
          {
            std::cout << GetName() << ": Bad step status combination for the same track " << std::endl;
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
        }
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
      //set the track ID
      m_Hit->set_trkid(aTrack->GetTrackID());
      m_SaveTrackId = aTrack->GetTrackID();
      //set the initial energy deposit
      m_Hit->set_edep(0);
      if (whichactive > 0)  // return of IsInOHCalDetector, > 0 hit in scintillator, < 0 hit in absorber
      {
        m_Hit->set_sector(sector_id);   // the sector id
        m_Hit->set_scint_id(tower_id);  // the slat id
        m_Hit->set_eion(0);
        m_Hit->set_raw_light_yield(0);  //  for scintillator only, initialize light yields
        m_Hit->set_light_yield(0);      //  for scintillator only, initialize light yields
        // Now save the container we want to add this hit to
        m_SaveHitContainer = m_HitContainer;
      }
      else
      {
        m_SaveHitContainer = m_AbsorberHitContainer;
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
      gSystem->Exit(1);
    }
    m_SavePostStepStatus = postPoint->GetStepStatus();
    // check if track id matches the initial one when the hit was created
    if (aTrack->GetTrackID() != m_SaveTrackId)
    {
      std::cout << GetName() << ": hits do not belong to the same track" << std::endl;
      std::cout << "saved track: " << m_SaveTrackId
                << ", current trackid: " << aTrack->GetTrackID()
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

    if (whichactive > 0)
    {
      m_Hit->set_eion(m_Hit->get_eion() + eion);
      light_yield = eion;
      if (m_LightScintModelFlag)
      {
        light_yield = GetVisibleEnergyDeposition(aStep);
        m_Hit->set_raw_light_yield(m_Hit->get_raw_light_yield() + light_yield);  // save raw Birks light yield
        if ((m_MapCorrHistChim[tower_id]) || (m_MapCorrHist[tower_id]))
        {
          const G4TouchableHandle& theTouchable = prePoint->GetTouchableHandle();
          const G4ThreeVector& worldPosition = postPoint->GetPosition();
          G4ThreeVector localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
          float lx = (localPosition.x() / cm);
          float ly = (localPosition.y() / cm);

          // convert to the map bin coordinates:
          int lcx = (int) (2.0 * lx) + 1;
          int lcy = (int) (2.0 * (ly + 0.5)) + 1;

          if ((sector_id == 29) || (sector_id == 30) || (sector_id == 31))
          {
            if ((lcy >= 1) && (lcy <= m_MapCorrHistChim[tower_id]->GetNbinsY()) &&
                (lcx >= 1) && (lcx <= m_MapCorrHistChim[tower_id]->GetNbinsX()))
            {
              light_yield *= (double) (m_MapCorrHistChim[tower_id]->GetBinContent(lcx, lcy));
            }
            else
            {
              light_yield = 0.0;
            }
          }
          else
          {
            if ((lcy >= 1) && (lcy <= m_MapCorrHist[tower_id]->GetNbinsY()) &&
                (lcx >= 1) && (lcx <= m_MapCorrHist[tower_id]->GetNbinsX()))
            {
              light_yield *= (double) (m_MapCorrHist[tower_id]->GetBinContent(lcx, lcy));
            }
            else
            {
              light_yield = 0.0;
            }
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
void PHG4OHCalSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_HitNodeName);
  m_AbsorberHitContainer = findNode::getClass<PHG4HitContainer>(topNode, m_AbsorberNodeName);

  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4OHCalSteppingAction::SetTopNode - unable to find " << m_HitNodeName << std::endl;
  }
  if (!m_AbsorberHitContainer)
  {
    if (Verbosity() > 1)
    {
      std::cout << "PHG4OHcalSteppingAction::SetTopNode - unable to find " << m_AbsorberNodeName << std::endl;
    }
  }
}

void PHG4OHCalSteppingAction::FieldChecker(const G4Step* aStep)
{
  Fun4AllServer* se = Fun4AllServer::instance();
  assert(se);

  static const std::string h_field_name = "hOHCalField";

  if (not se->isHistoRegistered(h_field_name))
  {
    TH2F* h = new TH2F(h_field_name.c_str(), "Magnetic field (Tesla) in HCal;X (cm);Y (cm)", 2400,
                       -300, 300, 2400, -300, 300);

    se->registerHisto(h, 1);

    std::cout << "PHG4OHCalSteppingAction::FieldChecker - make a histograme to check outer Hcal field map."
              << " Saved to Fun4AllServer Histo with name " << h_field_name << std::endl;
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

    const double Bz = FieldValue.z() / tesla;

    h->SetBinContent(binx, biny, Bz);

    std::cout << "PHG4OHCalSteppingAction::FieldChecker - "
              << "volume " << volume->GetName() << " / " << volume->GetLogicalVolume()->GetName()
              << "\t bin " << binx
              << ", " << biny << " : Bz= " << Bz << " B = " << FieldValue.mag() / tesla
              << " Tesla @ x,y = " << globPosVec[0] / cm
              << "," << globPosVec[1] / cm << " cm" << std::endl;
  }
}

void PHG4OHCalSteppingAction::SetHitNodeName(const std::string& type, const std::string& name)
{
  if (type == "G4HIT")
  {
    m_HitNodeName = name;
    return;
  }
  else if (type == "G4HIT_ABSORBER")
  {
    m_AbsorberNodeName = name;
    return;
  }
  std::cout << "Invalid output hit node type " << type << std::endl;
  gSystem->Exit(1);
  return;
}
