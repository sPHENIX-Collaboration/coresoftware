#include "PHG4MvtxSteppingAction.h"

#include "PHG4MvtxDetector.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Geant4/G4NavigationHistory.hh>
#include <Geant4/G4ParticleDefinition.hh>      // for G4ParticleDefinition
#include <Geant4/G4ReferenceCountedHandle.hh>  // for G4ReferenceCountedHandle
#include <Geant4/G4Step.hh>
#include <Geant4/G4StepPoint.hh>   // for G4StepPoint
#include <Geant4/G4StepStatus.hh>  // for fGeomBoundary, fAtRest...
#include <Geant4/G4String.hh>      // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>            // for G4ThreeVector
#include <Geant4/G4TouchableHandle.hh>        // for G4TouchableHandle
#include <Geant4/G4Track.hh>                  // for G4Track
#include <Geant4/G4TrackStatus.hh>            // for fStopAndKill
#include <Geant4/G4Types.hh>                  // for G4double
#include <Geant4/G4VPhysicalVolume.hh>        // for G4VPhysicalVolume
#include <Geant4/G4VTouchable.hh>             // for G4VTouchable
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

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

#include <cstdlib>  // for exit
#include <iostream>
#include <string>  // for operator<<, basic_string

class PHCompositeNode;

using namespace std;
//____________________________________________________________________________..
PHG4MvtxSteppingAction::PHG4MvtxSteppingAction(PHG4MvtxDetector* detector)
  : PHG4SteppingAction(detector->GetName())
  , m_Detector(detector)
  , m_HitContainer(nullptr)
  , m_AbsorberhitContainer(nullptr)
  , m_Hit(nullptr)
  , m_SaveShower(nullptr)
{
}

//____________________________________________________________________________..
PHG4MvtxSteppingAction::~PHG4MvtxSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete m_Hit;
}

//____________________________________________________________________________..
bool PHG4MvtxSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* sensor_volume = touch->GetVolume();

  // PHG4MvtxDetector->IsInMvtx(volume)
  // returns
  //  0 if outside of Mvtx
  //  1 if inside sensor

  // This checks if the volume is a sensor (doesn't tell us unique layer)
  // PHG4MvtxDetector->IsSensor(volume)
  // returns
  //  1 if volume is a sensor
  //  0 if not
  int whichactive = m_Detector->IsSensor(sensor_volume);

  if (!whichactive)
  {
    return false;
  }

  // This tells us if the volume belongs to the right stave for this layer
  // From the GDML file the 3rd volume up should be the half-stave
  // PHG4MvtxTelescopeDetector_->IsInMvtx(volume)
  // returns
  //  1 if in ladder belonging to this layer
  //  0 if not
  int layer_id = -9999;
  int stave_id = -9999;
  //cout << endl << "  In UserSteppingAction for layer " << layer_id << endl;
  G4VPhysicalVolume* vstave = touch->GetVolume(3);
  whichactive = m_Detector->IsInMvtx(vstave, layer_id, stave_id);
  if (layer_id < 0 || stave_id < 0)
  {
    cout << PHWHERE << "invalid Mvtx's layer (" << layer_id << ") or stave (" << stave_id << ") index " << endl;
    exit(1);
  }

  if (!whichactive)
  {
    return false;
  }

  if (Verbosity() > 5)
  {
    // make sure we know where we are!
    G4VPhysicalVolume* vtest = touch->GetVolume();
    cout << "Entering PHG4MvtxSteppingAction::UserSteppingAction for volume " << vtest->GetName() << endl;
    G4VPhysicalVolume* vtest1 = touch->GetVolume(1);
    cout << "Entering PHG4MvtxSteppingAction::UserSteppingAction for volume 1 up " << vtest1->GetName() << endl;
    G4VPhysicalVolume* vtest2 = touch->GetVolume(2);
    cout << "Entering PHG4MvtxSteppingAction::UserSteppingAction for volume 2 up " << vtest2->GetName() << endl;
    G4VPhysicalVolume* vtest3 = touch->GetVolume(3);
    cout << "Entering PHG4MvtxSteppingAction::UserSteppingAction for volume 3 up " << vtest3->GetName() << endl;
    G4VPhysicalVolume* vtest4 = touch->GetVolume(4);
    cout << "Entering PHG4MvtxSteppingAction::UserSteppingAction for volume 4 up " << vtest4->GetName() << endl;
  }

  //=======================================================================
  // We want the location of the hit
  // Here we will record in the hit object:
  //   The stave, half-stave, module and chip numbers
  //   The energy deposited
  //   The entry point and exit point in world coordinates
  //   The entry point and exit point in local (sensor) coordinates
  // The pixel number will be derived later from the entry and exit points in the sensor local coordinates
  //=======================================================================

  //int stave_number = -1; // stave_number from stave map values
  int half_stave_number = -1;
  int module_number = -1;
  int chip_number = -1;

  if (Verbosity() > 0)
    cout << endl
         << "  UserSteppingAction: layer " << layer_id;
  boost::char_separator<char> sep("_");
  boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;

  //OLD ITS.gdml: chip number is from  "ITSUChip[layer number]_[chip number]
  //NEW: chip number is from  "MVTXChip_[chip number]
  G4VPhysicalVolume* v1 = touch->GetVolume(1);
  boost::tokenizer<boost::char_separator<char> > tok1(v1->GetName(), sep);
  tokeniter = tok1.begin();
  ++tokeniter;
  chip_number = boost::lexical_cast<int>(*tokeniter);
  if (Verbosity() > 0) cout << " chip  " << chip_number;
  G4VPhysicalVolume* v2 = touch->GetVolume(2);
  boost::tokenizer<boost::char_separator<char> > tok2(v2->GetName(), sep);
  tokeniter = tok2.begin();
  ++tokeniter;
  module_number = boost::lexical_cast<int>(*tokeniter);
  if (Verbosity() > 0) cout << " module " << module_number;

  // The stave number  is the imprint number from the assembly volume imprint
  // The assembly volume history string format is (e.g.):
  // av_13_impr_4_ITSUHalfStave6_pv_1
  // where "impr_4" says stave number 4
  // and where "ITSUHalfStave6_pv_1" says hald stave number 1 in layer number 6

  G4VPhysicalVolume* v3 = touch->GetVolume(3);
  boost::tokenizer<boost::char_separator<char> > tok3(v3->GetName(), sep);
  tokeniter = tok3.begin();
  ++tokeniter;
  ++tokeniter;
  ++tokeniter;
  //stave_number = boost::lexical_cast<int>(*tokeniter) - 1;  // starts counting imprints at 1, we count staves from 0!
  if (Verbosity() > 0) cout << " stave " << stave_id;
  ++tokeniter;
  ++tokeniter;
  ++tokeniter;
  half_stave_number = boost::lexical_cast<int>(*tokeniter);
  if (Verbosity() > 0) cout << " half_stave " << half_stave_number;

  // FYI: doing string compares inside a stepping action sounds like a recipe
  // for failure inside a heavy ion event... we'll wait and see how badly
  // this profiles. -MPM

  // Now we want to collect information about the hit

  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  const G4Track* aTrack = aStep->GetTrack();
  if (Verbosity() > 0) cout << " edep = " << edep << endl;

  // if this cylinder stops everything, just put all kinetic energy into edep
  if (m_Detector->IsBlackHole(layer_id))
  {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track*>(aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  // test if we are active
  if (m_Detector->IsActive(layer_id))
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

    G4ThreeVector worldPosition;  // localPosition;
    G4TouchableHandle theTouchable;

    G4StepPoint* prePoint = aStep->GetPreStepPoint();
    G4StepPoint* postPoint = aStep->GetPostStepPoint();

    if (Verbosity() > 0)
    {
      G4ParticleDefinition* def = aTrack->GetDefinition();
      cout << " Particle: " << def->GetParticleName() << endl;
    }

    switch (prePoint->GetStepStatus())
    {
    case fGeomBoundary:
    case fUndefined:
      if (!m_Hit)
      {
        m_Hit = new PHG4Hitv1();
      }
      m_Hit->set_layer((unsigned int) layer_id);

      // set the index values needed to locate the sensor strip
      m_Hit->set_property(PHG4Hit::prop_stave_index, stave_id);
      m_Hit->set_property(PHG4Hit::prop_half_stave_index, half_stave_number);
      m_Hit->set_property(PHG4Hit::prop_module_index, module_number);
      m_Hit->set_property(PHG4Hit::prop_chip_index, chip_number);

      worldPosition = prePoint->GetPosition();

      if (Verbosity() > 0)
      {
        theTouchable = prePoint->GetTouchableHandle();
        cout << "entering: depth = " << theTouchable->GetHistory()->GetDepth() << endl;
        G4VPhysicalVolume* vol1 = theTouchable->GetVolume();
        cout << "entering volume name = " << vol1->GetName() << endl;
      }

      /*
	  localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
	  */

      // Store the local coordinates for the entry point
      StoreLocalCoordinate(m_Hit, aStep, true, false);

      // Store the entrance values in cm in world coordinates
      m_Hit->set_x(0, prePoint->GetPosition().x() / cm);
      m_Hit->set_y(0, prePoint->GetPosition().y() / cm);
      m_Hit->set_z(0, prePoint->GetPosition().z() / cm);

      m_Hit->set_px(0, prePoint->GetMomentum().x() / GeV);
      m_Hit->set_py(0, prePoint->GetMomentum().y() / GeV);
      m_Hit->set_pz(0, prePoint->GetMomentum().z() / GeV);

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

      // Now add the hit
      //	  m_Hit->print();
      //m_HitContainer->AddHit(layer_id, m_Hit);
      break;
    default:
      break;
    }

    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist

    /*
      // Note that the after you reach the boundary the touchable for the postPoint points to the next volume, not the sensor!
      // This was given to me as the way to get back to the sensor volume, but it does not work
      theTouchable = postPoint->GetTouchableHandle();
      localPosition = history->GetTransform(history->GetDepth() - 1).TransformPoint(worldPosition);
      cout << "Exit local coords: x " <<  localPosition.x() / cm << " y " <<  localPosition.y() / cm << " z " <<  localPosition.z() / cm << endl;
      */

    /*
      // Use the prePoint from the final step  for now, until I understand how to get the exit point in the sensor volume coordinates
      //======================================================================================
      theTouchable = prePoint->GetTouchableHandle();
      vol2 = theTouchable->GetVolume();
      if(Verbosity() > 0)
	cout << "exiting volume name = " << vol2->GetName() << endl;
      worldPosition = prePoint->GetPosition();
      if(Verbosity() > 0)
	cout << "Exit world coords prePoint: x " <<  worldPosition.x() / cm << " y " <<  worldPosition.y() / cm << " z " <<  worldPosition.z() / cm << endl;

      // This is for consistency with the local coord position, the world coordinate exit position is correct
      m_Hit->set_x( 1, prePoint->GetPosition().x() / cm );
      m_Hit->set_y( 1, prePoint->GetPosition().y() / cm );
      m_Hit->set_z( 1, prePoint->GetPosition().z() / cm );

      const G4NavigationHistory *history = theTouchable->GetHistory();
      //cout << "exiting: depth = " << history->GetDepth() <<  " volume name = " << history->GetVolume(history->GetDepth())->GetName() << endl;
      localPosition = history->GetTransform(history->GetDepth()).TransformPoint(worldPosition);

      m_Hit->set_local_x(1, localPosition.x() / cm);
      m_Hit->set_local_y(1, localPosition.y() / cm);
      m_Hit->set_local_z(1, localPosition.z() / cm);
      */

    // Store the local coordinates for the exit point
    StoreLocalCoordinate(m_Hit, aStep, false, true);

    // Store world coordinates for the exit point
    m_Hit->set_x(1, postPoint->GetPosition().x() / cm);
    m_Hit->set_y(1, postPoint->GetPosition().y() / cm);
    m_Hit->set_z(1, postPoint->GetPosition().z() / cm);

    m_Hit->set_px(1, postPoint->GetMomentum().x() / GeV);
    m_Hit->set_py(1, postPoint->GetMomentum().y() / GeV);
    m_Hit->set_pz(1, postPoint->GetMomentum().z() / GeV);

    m_Hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
    //sum up the energy to get total deposited
    m_Hit->set_edep(m_Hit->get_edep() + edep);
    if (geantino)
    {
      m_Hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
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

    if (Verbosity() > 0)
    {
      G4StepPoint* prePointA = aStep->GetPreStepPoint();
      G4StepPoint* postPointA = aStep->GetPostStepPoint();
      cout << "----- PHg4MvtxSteppingAction::UserSteppingAction - active volume = " << sensor_volume->GetName() << endl;
      cout << "       layer = " << layer_id << endl;
      cout << "       stave number = " << stave_id << " half_stave_number = " << half_stave_number << endl;
      cout << "       module number  = " << module_number << endl;
      cout << "       chip number = " << chip_number << endl;
      cout << "       prepoint x position " << prePointA->GetPosition().x() / cm << endl;
      cout << "       prepoint y position " << prePointA->GetPosition().y() / cm << endl;
      cout << "       prepoint z position " << prePointA->GetPosition().z() / cm << endl;
      cout << "       postpoint x position " << postPointA->GetPosition().x() / cm << endl;
      cout << "       postpoint y position " << postPointA->GetPosition().y() / cm << endl;
      cout << "       postpoint z position " << postPointA->GetPosition().z() / cm << endl;
      cout << "       edep " << edep << endl;
    }

    if (Verbosity() > 0)
    {
      cout << "  stepping action found hit:" << endl;
      m_Hit->print();
      cout << endl
           << endl;
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
        m_HitContainer->AddHit(layer_id, m_Hit);
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

    // return true to indicate the hit was used
    return true;
  }
  else
  {
    return false;
  }

  return false;
}

//____________________________________________________________________________..
void PHG4MvtxSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  m_HitContainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  m_AbsorberhitContainer = findNode::getClass<PHG4HitContainer>(topNode, absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!m_HitContainer)
  {
    std::cout << "PHG4MvtxSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if (!m_AbsorberhitContainer)
  {
    if (Verbosity() > 0)
    {
      cout << "PHG4MvtxSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}
