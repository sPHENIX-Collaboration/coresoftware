#include "PHG4SiliconTrackerSteppingAction.h"
#include "PHG4SiliconTrackerDefs.h"
#include "PHG4SiliconTrackerDetector.h"

#include "PHG4StepStatusDecode.h"

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

#include <boost/foreach.hpp>
#include <boost/format.hpp>
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

  // set ladder index
  int sphxlayer = 0;
  int inttlayer = 0;
  int ladderz = 0;
  int ladderphi = 0;
  int zposneg = 0;
  int strip_z_index = 0;
  int strip_y_index = 0;
  // int save_strip_z_index = 0;
  // int save_strip_y_index = 0;
  // int save_ladderz = 0;
  // int save_ladderphi = 0;
  if (whichactive > 0)  // silicon active sensor
  {
    if (Verbosity() > 0)
    {
      cout << endl
           << "PHG4SilicoTrackerSteppingAction::UserSteppingAction for volume name (pre) " << touch->GetVolume()->GetName()
           << " volume name (2) " << touch->GetVolume(2)->GetName()
           << " volume->GetTranslation " << touch->GetVolume()->GetTranslation()
           << " volume->GetCopyNo() " << volume->GetCopyNo()
           << endl;
    }

    // Get the layer and ladder information which are 2 steps up in the volume hierarchy
    // the ladder also contains inactive volumes but we check in m_Detector->IsInSiliconTracker(volume)
    // if we are in an active logical volume whioch is located in this ladder
    auto iter = m_Detector->get_ActiveVolumeTuple(touch->GetVolume(2));
    tie(inttlayer, ladderz, ladderphi, zposneg) = iter->second;
    if (inttlayer < 0 || inttlayer > 3)
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

    // Find the strip y and z index values from the copy number (integer division, quotient is strip_y, remainder is strip_z)
    int laddertype = (m_LadderTypeMap.find(inttlayer))->second;
    int nstrips_z_sensor;
    switch (ladderz)
    {
    case 0:
      nstrips_z_sensor = m_nStripsZSensor.find(laddertype)->second.first;
      break;
    case 1:
      nstrips_z_sensor = m_nStripsZSensor.find(laddertype)->second.second;
      break;
    default:
      cout << "Bad ladderz: " << ladderz << endl;
      gSystem->Exit(1);
      // this is just to make the optimizer happy which otherwise complains about possibly
      // uninitialized variables. It doesn't know gSystem->Exit(1) quits,
      // this exit here terminates the program for it
      exit(1);
    }
    div_t copydiv = div(volume->GetCopyNo(), nstrips_z_sensor);
    strip_y_index = copydiv.quot;
    strip_z_index = copydiv.rem;
    // save_strip_z_index = strip_z_index;
    // save_strip_y_index = strip_y_index;
    // save_ladderz = ladderz;
    // save_ladderphi = ladderphi;
    G4ThreeVector strip_pos = volume->GetTranslation();
    G4ThreeVector prepos = prePoint->GetPosition();
    G4ThreeVector postpos = postPoint->GetPosition();
    if (Verbosity() > 0)
    {
      cout << " sphxlayer " << sphxlayer << " inttlayer " << inttlayer << " ladderz " << ladderz << " ladderphi " << ladderphi
           << " zposneg " << zposneg
           << " copy no. " << volume->GetCopyNo() << " nstrips_z_sensor " << nstrips_z_sensor
           << " strip_y_index " << strip_y_index << " strip_z_index " << strip_z_index << endl;
    }
    // There are two failure modes observed for this stupid parameterised volume:
    //  1) If the prePoint step status is "fUndefined" then the copy number is sometimes kept from the last hit, which is often an unrelated volume
    //  2) If the  pre and post step are in the same volume but they both have status fGeomBoundary, the volume is assigned the copy number of the next volume, which is usually off by one in strip_y_index
    // in both cases we need to find the correct strip_y_index and strip_z_index values the hard way - that is what is done here
    int fixit = 0;
    G4VPhysicalVolume* volume_post = postPoint->GetTouchableHandle()->GetVolume();
    //      G4VPhysicalVolume* volume_pre = prePoint->GetTouchableHandle()->GetVolume();
    G4LogicalVolume* logvolpre = volume->GetLogicalVolume();
    G4LogicalVolume* logvolpost = volume_post->GetLogicalVolume();
    if (prePoint->GetStepStatus() == fGeomBoundary && postPoint->GetStepStatus() == fGeomBoundary)
    {
      // this is just failsafe - we still have those impossible hits where pre and poststep
      // are in the same volume with status fGeomBoundary
      // but the extraction of the strip index above works for those
      // my assumption is that the end of a physics step coincides with the boundary and
      // the physics step size is taken rather the geometric step size
      if (logvolpre == logvolpost)
      {
        if (volume->GetCopyNo() == volume_post->GetCopyNo())
        {
          fixit = 1;
        }
      }
    }

    if (prePoint->GetStepStatus() == fUndefined)
    {
      fixit = 2;

      // cout << "copy pre vol: " << volume->GetCopyNo() << ", " << volume->GetName()
      // 	   << ", post: " << volume_post->GetCopyNo() << ", " << volume_post->GetName()
      // 	   << ", pre: " << volume_pre->GetCopyNo() << ", " << volume_pre->GetName()
      // 	   << endl;
      // if (prePoint->GetProcessDefinedStep())
      // {
      // cout << prePoint->GetProcessDefinedStep()->GetProcessName() << endl;
      // }
      // else
      // {
      // 	cout << "Process is null ptr, try post step" << endl;
      // 	if (postPoint->GetProcessDefinedStep())
      // 	{
      // 	  cout << postPoint->GetProcessDefinedStep()->GetProcessName() << endl;
      // 	}
      // 	else
      // 	{
      // 	  cout << "no post step proc either" << endl;
      // 	}
      // }
      // G4ThreeVector preworldPos = prePoint->GetPosition();
      // G4ThreeVector strip_pos = touch->GetHistory()->GetTransform(touch->GetHistory()->GetDepth() - 1).TransformPoint(preworldPos);
      // G4ThreeVector postworldPos = postPoint->GetPosition();
      // G4ThreeVector poststrip_pos = touch->GetHistory()->GetTransform(touch->GetHistory()->GetDepth() - 1).TransformPoint(postworldPos);
      // cout << "preworldpos: " << preworldPos << endl;
      // cout << "pre strip pos: " << strip_pos << endl;
      // cout << "postworldPos: " << postPoint->GetPosition() << endl;
      // cout << "poststrip_pos: " << poststrip_pos << endl;
    }

    if (fixit)
    {
      //	cout << "fixit: " << fixit << endl;
      int strip_y_index_old = strip_y_index;
      int strip_z_index_old = strip_z_index;

      // we need a hack to compare the values above with the correct strip index values
      // the transform of the world coordinates into the sensor frame will work correctly,
      // so we determine the strip indices from the hit position
      G4ThreeVector preworldPos = prePoint->GetPosition();
      G4ThreeVector strip_pos = touch->GetHistory()->GetTransform(touch->GetHistory()->GetDepth() - 1).TransformPoint(preworldPos);
      G4ThreeVector postworldPos = postPoint->GetPosition();
      G4ThreeVector poststrip_pos = touch->GetHistory()->GetTransform(touch->GetHistory()->GetDepth() - 1).TransformPoint(postworldPos);

      strip_z_index = 0;
      double strip_z;
      int nstrips_z_sensor;
      switch (ladderz)
      {
      case 0:
        strip_z = m_StripZMap.find(laddertype)->second.first;
        nstrips_z_sensor = m_nStripsZSensor.find(laddertype)->second.first;
        break;
      case 1:
        strip_z = m_StripZMap.find(laddertype)->second.second;
        nstrips_z_sensor = m_nStripsZSensor.find(laddertype)->second.second;
        break;
      default:
        cout << "Bad ladderz: " << ladderz << endl;
        gSystem->Exit(1);
        // this is just to make the optimizer happy which otherwise complains about possibly
        // uninitialized variables. It doesn't know gSystem->Exit(1) quits,
        // this exit here terminates the program for it
        exit(1);
      }
      // cout << "nstrips_z_sensor: " << nstrips_z_sensor
      //      << ", strip_z: " << strip_z
      //      << endl;
      for (int i = 0; i < nstrips_z_sensor; ++i)
      {
        const double zmin = strip_z * i - (strip_z / 2.) * nstrips_z_sensor;
        const double zmax = strip_z * (i + 1) - (strip_z / 2.) * nstrips_z_sensor;
        if (strip_pos.z() / mm > zmin && strip_pos.z() / mm <= zmax)
        {
          strip_z_index = i;
          if (Verbosity() > 1) std::cout << "                            revised strip z position = " << strip_z_index << std::endl;
          break;
        }
      }

      strip_y_index = 0;
      int nstrips_phi_cell = m_nStripsPhiCell.find(laddertype)->second;
      double strip_y = m_StripYMap.find(laddertype)->second;
      // cout << "nstrips_phi_cell: " << nstrips_phi_cell
      //      << ", strip_y: " << strip_y
      //      << endl;
      for (int i = 0; i < 2 * nstrips_phi_cell; ++i)
      {
        const double ymin = strip_y * i - strip_y * nstrips_phi_cell;
        ;
        const double ymax = strip_y * (i + 1) - strip_y * nstrips_phi_cell;
        if (strip_pos.y() / mm > ymin && strip_pos.y() / mm <= ymax)
        {
          strip_y_index = i;
          if (Verbosity() > 1) std::cout << "                            revised strip y position = " << strip_y_index << std::endl;
          break;
        }
      }

      if (fixit == 1)
      {
        if (strip_y_index_old != strip_y_index || strip_z_index_old != strip_z_index)
        {
          G4VPhysicalVolume* volume_post = postPoint->GetTouchableHandle()->GetVolume();
          G4LogicalVolume* logvolpre = volume->GetLogicalVolume();
          G4LogicalVolume* logvolpost = volume_post->GetLogicalVolume();
          G4ThreeVector preworldPos = prePoint->GetPosition();
          G4ThreeVector strip_pos = touch->GetHistory()->GetTransform(touch->GetHistory()->GetDepth() - 1).TransformPoint(preworldPos);
          G4ThreeVector postworldPos = postPoint->GetPosition();
          G4ThreeVector poststrip_pos = touch->GetHistory()->GetTransform(touch->GetHistory()->GetDepth() - 1).TransformPoint(postworldPos);

          cout << "Overlap detected in volume " << volume->GetName() << " where post volume "
               << volume_post->GetName() << " has same copy no." << volume->GetCopyNo()
               << " pre and post step point of same volume for step status fGeomBoundary" << endl;
          cout << "logvol name " << logvolpre->GetName() << ", post: " << logvolpost->GetName() << endl;
          cout << "strip y bef: " << strip_y_index_old << ", strip z: " << strip_z_index_old << endl;
          cout << " strip y aft: " << strip_y_index << ", strip z: " << strip_z_index << endl;
          cout << "pre hitpos x: " << strip_pos.x() << ", y: " << strip_pos.y() << ", z: "
               << strip_pos.z() << endl;
          cout << "posthitpos x: " << poststrip_pos.x() << ", y: " << poststrip_pos.y() << ", z: "
               << poststrip_pos.z() << endl;
          cout << "eloss: " << aStep->GetTotalEnergyDeposit() / GeV << " GeV" << endl;
          cout << "safety prestep: " << prePoint->GetSafety()
               << ", poststep: " << postPoint->GetSafety() << endl;
        }
      }
      else
      {
        if (Verbosity() > 1) cout << "Detected fUndefined for prePoint step status, re-calculated strip_z_index and strip_y_index independently above" << endl;
      }
    }
  }     // end of whichactive > 0 block
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
    cout << "prePoint step status = " << prePoint->GetStepStatus() << " postPoint step status = " << postPoint->GetStepStatus() << endl;
  }
  switch (prePoint->GetStepStatus())
  {
  case fGeomBoundary:
  case fUndefined:

    if (Verbosity() > 1)
    {
      cout << " found prePoint step status of fGeomBoundary or fUndefined, start a new hit " << endl;
    }
    // if previous hit was saved, hit pointer was set to nullptr
    // and we have to make a new one
    if (!m_Hit)
    {
      m_Hit = new PHG4Hitv1();
    }

    m_Hit->set_layer((unsigned int) sphxlayer);

    // set the index values needed to locate the sensor strip
    if (zposneg == 1) ladderz += 2;  // ladderz = 0, 1 for negative z and = 2, 3 for positive z
    m_Hit->set_ladder_z_index(ladderz);

    if (whichactive > 0)
    {
      m_Hit->set_strip_z_index(strip_z_index);
      m_Hit->set_strip_y_index(strip_y_index);
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

    /*
    cout << "     hit position x,y,z = " << prePoint->GetPosition().x() / cm
    << "    " << prePoint->GetPosition().y() / cm
    << "     " << prePoint->GetPosition().z() / cm
    << endl;
    */

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

  // here we just update the exit values, it will be overwritten
  // for every step until we leave the volume or the particle
  // ceases to exist
  m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
  m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
  m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

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
      cout << " postPoint step status changed to " << postPoint->GetStepStatus() << " save hit and delete it" << endl;
      cout << " fWorldBoundary " << fWorldBoundary
           << " fGeomBoundary = " << fGeomBoundary
           << " fAtRestDoItProc " << fAtRestDoItProc
           << " fAlongStepDoItProc " << fAlongStepDoItProc
           << " fPostStepDoItProc " << fPostStepDoItProc
           << " fUserDefinedLimit " << fUserDefinedLimit
           << " fExclusivelyForcedProc " << fExclusivelyForcedProc
           << " fUndefined " << fUndefined
           << " fStopAndKill " << fStopAndKill
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

  if (Verbosity() > 1)
  {
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();
    G4ThreeVector preworldPos = prePoint->GetPosition();
    G4ThreeVector postworldPos = postPoint->GetPosition();

    cout << " entry point world pos " << prePoint->GetPosition().x() << "  " << prePoint->GetPosition().y() << "  " << prePoint->GetPosition().z() << endl;
    cout << " exit point world pos " << postPoint->GetPosition().x() << "  " << postPoint->GetPosition().y() << "  " << postPoint->GetPosition().z() << endl;

    // The exit point transforms do not work here because the particle has already entered the next volume
    // - go back and find Jin's fix for this

    // strip local pos
    G4TouchableHistory* theTouchable = (G4TouchableHistory*) (prePoint->GetTouchable());
    G4ThreeVector prelocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(preworldPos);
    cout << " entry point strip local pos: "
         << " is " << prelocalPos.x() << " " << prelocalPos.y() << " " << prelocalPos.z() << endl;
    G4TouchableHistory* postTouchable = (G4TouchableHistory*) (postPoint->GetTouchable());
    G4ThreeVector postlocalPos = postTouchable->GetHistory()->GetTopTransform().TransformPoint(postworldPos);
    cout << " exit point strip local pos: " << postlocalPos.x() << " " << postlocalPos.y() << " " << postlocalPos.z() << endl;

    // sensor local pos
    G4ThreeVector presensorLocalPos = theTouchable->GetHistory()->GetTransform(theTouchable->GetHistory()->GetDepth() - 1).TransformPoint(preworldPos);
    cout << " entry point sensor local pos: " << presensorLocalPos.x() << " " << presensorLocalPos.y() << " " << presensorLocalPos.z() << endl;
    G4ThreeVector postsensorLocalPos = postTouchable->GetHistory()->GetTransform(postTouchable->GetHistory()->GetDepth() - 1).TransformPoint(postworldPos);
    cout << " exit point sensor local pos: " << postsensorLocalPos.x() << " " << postsensorLocalPos.y() << " " << postsensorLocalPos.z() << endl;
  }

  if (whichactive > 0 && Verbosity() > 0)  // return of IsInSiliconTracker, > 0 hit in si-strip, < 0 hit in absorber
  {
    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();
    cout << "----- PHg4SiliconTrackerSteppingAction::UserSteppingAction - active volume = " << volume->GetName() << endl;
    cout << "       strip_z_index = " << strip_z_index << " strip_y_index = " << strip_y_index << endl;
    cout << "       prepoint x position " << prePoint->GetPosition().x() / cm << endl;
    cout << "       prepoint y position " << prePoint->GetPosition().y() / cm << endl;
    cout << "       prepoint z position " << prePoint->GetPosition().z() / cm << endl;
    cout << "       postpoint x position " << postPoint->GetPosition().x() / cm << endl;
    cout << "       postpoint y position " << postPoint->GetPosition().y() / cm << endl;
    cout << "       postpoint z position " << postPoint->GetPosition().z() / cm << endl;
    cout << "       edep " << edep << endl;
    cout << "       eion " << eion << endl;
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
