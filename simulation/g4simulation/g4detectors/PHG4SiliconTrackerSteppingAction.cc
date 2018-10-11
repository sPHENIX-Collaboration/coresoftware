#include "PHG4SiliconTrackerSteppingAction.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4SiliconTrackerDefs.h"
#include "PHG4SiliconTrackerDetector.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <TSystem.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4TouchableHistory.hh>
#include <Geant4/G4VProcess.hh>

#include <gsl/gsl_math.h>

#include <iostream>

using namespace std;

//____________________________________________________________________________..
PHG4SiliconTrackerSteppingAction::PHG4SiliconTrackerSteppingAction(PHG4SiliconTrackerDetector* detector, const PHParametersContainer* parameters, const pair<vector<pair<int, int>>::const_iterator, vector<pair<int, int>>::const_iterator>& layer_begin_end)
  : m_Detector(detector)
  , m_Hits(nullptr)
  , m_AbsorberHits(nullptr)
  , m_Hit(nullptr)
  , m_SaveHitContainer(nullptr)
  , m_SaveShower(nullptr)
  , m_ParamsContainer(parameters)
{
  // loop over layers to get laddertype nd active status for each layer
  for (auto layeriter = layer_begin_end.first; layeriter != layer_begin_end.second; ++layeriter)
  {
    int layer = layeriter->second;
    const PHParameters* par = m_ParamsContainer->GetParameters(layer);
    m_IsActiveMap[layer] = par->get_int_param("active");
    m_IsBlackHoleMap[layer] = par->get_int_param("blackhole");
    m_LadderTypeMap.insert(make_pair(layer, par->get_int_param("laddertype")));
    m_InttToTrackerLayerMap.insert(make_pair(layeriter->second, layeriter->first));
  }

  // Get the parameters for each laddertype
  for (auto iter = PHG4SiliconTrackerDefs::m_SensorSegmentationSet.begin(); iter != PHG4SiliconTrackerDefs::m_SensorSegmentationSet.end(); ++iter)
  {
    const PHParameters* par = m_ParamsContainer->GetParameters(*iter);
    m_StripYMap.insert(make_pair(*iter, par->get_double_param("strip_y") * cm));
    m_StripZMap.insert(make_pair(*iter, make_pair(par->get_double_param("strip_z_0") * cm, par->get_double_param("strip_z_1") * cm)));
    m_nStripsPhiCell.insert(make_pair(*iter, par->get_int_param("nstrips_phi_cell")));
    m_nStripsZSensor.insert(make_pair(*iter, make_pair(par->get_int_param("nstrips_z_sensor_0"), par->get_int_param("nstrips_z_sensor_1"))));
  }
}

PHG4SiliconTrackerSteppingAction::~PHG4SiliconTrackerSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4SiliconTrackerSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();
  const G4Track* aTrack = aStep->GetTrack();
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();

  const int whichactive = m_Detector->IsInSiliconTracker(volume);

  if (!whichactive)
  {
    return false;
  }

  if (Verbosity() > 0)
  {
    cout << endl
         << "PHG4SilicoTrackerSteppingAction::UserSteppingAction for volume name (pre) " << touch->GetVolume()->GetName()
         << " volume name (1) " << touch->GetVolume(1)->GetName()
         << " volume->GetTranslation " << touch->GetVolume()->GetTranslation()
         << endl;
  }

  // set ladder index
  int sphxlayer = 0;
  int inttlayer = 0;
  int ladderz = 0;
  int ladderphi = 0;
  int zposneg = 0;

  if (whichactive > 0)  // silicon active sensor
  {
    // Get the layer and ladder information which is step up in the volume hierarchy
    // the ladder also contains inactive volumes but we check in m_Detector->IsInSiliconTracker(volume)
    // if we are in an active logical volume whioch is located in this ladder
    auto iter = m_Detector->get_ActiveVolumeTuple(touch->GetVolume(1));
    tie(inttlayer, ladderz, ladderphi, zposneg) = iter->second;
    if (inttlayer < 0 || inttlayer > 7)
    {
      assert(!"PHG4SiliconTrackerSteppingAction: check INTT ladder layer.");
    }
    sphxlayer = m_InttToTrackerLayerMap.find(inttlayer)->second;
    map<int, int>::const_iterator activeiter = m_IsActiveMap.find(inttlayer);
    if (activeiter == m_IsActiveMap.end())
    {
      cout << "PHG4SiliconTrackerSteppingAction: could not find active flag for layer " << inttlayer << endl;
      gSystem->Exit(1);
    }
    if (activeiter->second == 0)
    {
      return false;
    }
  }
  else  // whichactive < 0, silicon inactive area, FPHX, stabe etc. as absorbers
  {
    auto iter = m_Detector->get_PassiveVolumeTuple(touch->GetVolume(0)->GetLogicalVolume());
    tie(inttlayer, ladderz) = iter->second;
    sphxlayer = inttlayer;  //for absorber we use the INTT layer, not the tracking layer in sPHENIX
  }                         // end of si inactive area block

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  // if this block stops everything, just put all kinetic energy into edep
  if ((m_IsBlackHoleMap.find(inttlayer))->second == 1)
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  bool geantino = false;

  // the check for the pdg code speeds things up, I do not want to make
  // an expensive string compare for every track when we know
  // geantino or chargedgeantino has pid=0
  if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 && aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
  {
    geantino = true;
  }

  if (Verbosity() > 1)
  {
    cout << " aTrack->GetTrackID " << aTrack->GetTrackID() << " aTrack->GetParentID " << aTrack->GetParentID()
         << " INTT layer " << inttlayer
         << " prePoint step status = " << prePoint->GetStepStatus() << " postPoint step status = " << postPoint->GetStepStatus()
         << endl;
  }
  switch (prePoint->GetStepStatus())
  {
  case fGeomBoundary:
  case fUndefined:

    if (Verbosity() > 1)
    {
      cout << "  start new g4hit for:  aTrack->GetParentID " << aTrack->GetParentID() << " aTrack->GetTrackID " << aTrack->GetTrackID() << " INTT layer " << inttlayer
           << "   prePoint step status  = " << prePoint->GetStepStatus() << " postPoint->GetStepStatus = " << postPoint->GetStepStatus() << endl;
    }
    // if previous hit was saved, hit pointer was set to nullptr
    // and we have to make a new one
    if (!m_Hit)
    {
      m_Hit = new PHG4Hitv1();
    }

    // set the index values needed to locate the sensor strip
    if (zposneg == 1)
    {
      ladderz += 2;  // ladderz = 0, 1 for negative z and = 2, 3 for positive z
    }

    m_Hit->set_ladder_z_index(ladderz);

    if (whichactive > 0)
    {
      m_Hit->set_layer(sphxlayer);
      m_Hit->set_ladder_phi_index(ladderphi);
      m_Hit->set_px(0, prePoint->GetMomentum().x() / GeV);
      m_Hit->set_py(0, prePoint->GetMomentum().y() / GeV);
      m_Hit->set_pz(0, prePoint->GetMomentum().z() / GeV);
      m_Hit->set_eion(0);
    }

    //here we set the entrance values in cm
    m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
    m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
    m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

    StoreLocalCoordinate(m_Hit, aStep, true, false);

    if (Verbosity() > 1)
    {
      cout << "     prePoint hit position x,y,z = " << prePoint->GetPosition().x() / cm
           << "    " << prePoint->GetPosition().y() / cm
           << "     " << prePoint->GetPosition().z() / cm
           << endl;
    }

    // time in ns
    m_Hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);

    //set the track ID
    m_Hit->set_trkid(aTrack->GetTrackID());

    //set the initial energy deposit
    m_Hit->set_edep(0);

    if (whichactive > 0)  // return of IsInSiliconTracker, > 0 hit in si-strip, < 0 hit in absorber
    {
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

  //cout << "  Update exit values for prePoint->GetStepStatus = " << prePoint->GetStepStatus() << " and postPoint->GetStepStatus = " << postPoint->GetStepStatus() << endl;

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

  if (geantino)
  {
    m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
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
  if (postPoint->GetStepStatus() == fGeomBoundary ||
      postPoint->GetStepStatus() == fWorldBoundary ||
      postPoint->GetStepStatus() == fAtRestDoItProc ||
      aTrack->GetTrackStatus() == fStopAndKill)
  {
    if (Verbosity() > 1)
    {
      cout << " postPoint step status changed to " << postPoint->GetStepStatus() << " or aTrack status changed to "
           << aTrack->GetTrackStatus() << endl;
      cout << "  end g4hit for:  aTrack->GetParentID " << aTrack->GetParentID() << " aTrack->GetTrackID " << aTrack->GetTrackID() << " eion = " << eion << endl;
      cout << "     end hit position x,y,z = " << postPoint->GetPosition().x() / cm
           << "    " << postPoint->GetPosition().y() / cm
           << "     " << postPoint->GetPosition().z() / cm
           << endl;
    }

    // save only hits with energy deposit (or -1 for geantino)
    if (m_Hit->get_edep())
    {
      m_SaveHitContainer->AddHit(sphxlayer, m_Hit);
      if (m_SaveShower)
      {
        m_SaveShower->add_g4hit_id(m_SaveHitContainer->GetID(), m_Hit->get_hit_id());
      }
      if (Verbosity() > 1)
      {
        m_Hit->print();
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

//____________________________________________________________________________..
void PHG4SiliconTrackerSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  const string detectorname = (m_Detector->SuperDetector() != "NONE") ? m_Detector->SuperDetector() : m_Detector->GetName();
  const string hitnodename = "G4HIT_" + detectorname;
  const string absorbernodename = "G4HIT_ABSORBER_" + detectorname;

  //now look for the map and grab a pointer to it.
  m_Hits = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  m_AbsorberHits = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!m_Hits)
    cout << "PHG4SiliconTrackerSteppingAction::SetTopNode - unable to find " << hitnodename << endl;

  if (!m_AbsorberHits && Verbosity() > 1)
    cout << "PHG4SiliconTrackerSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
}
