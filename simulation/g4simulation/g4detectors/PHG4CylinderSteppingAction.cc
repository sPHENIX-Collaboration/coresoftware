#include "PHG4CylinderSteppingAction.h"
#include "PHG4CylinderDetector.h"
#include "PHG4StepStatusDecode.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4UserLimits.hh>

#include <boost/io/ios_state.hpp>
#include <iomanip>
#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4CylinderSteppingAction::PHG4CylinderSteppingAction(PHG4CylinderDetector* detector, const PHParameters* parameters)
  : detector_(detector)
  , params(parameters)
  , hits_(nullptr)
  , hit(nullptr)
  , saveshower(nullptr)
  , savevolpre(nullptr)
  , savevolpost(nullptr)
  , save_light_yield(params->get_int_param("lightyield"))
  , savetrackid(-1)
  , saveprestepstatus(-1)
  , savepoststepstatus(-1)
  , active(params->get_int_param("active"))
  , IsBlackHole(params->get_int_param("blackhole"))
  , use_g4_steps(params->get_int_param("use_g4steps"))
  , zmin(params->get_double_param("place_z") * cm - params->get_double_param("length") * cm / 2.)
  , zmax(params->get_double_param("place_z") * cm + params->get_double_param("length") * cm / 2.)
  , tmin(params->get_double_param("tmin") * ns)
  , tmax(params->get_double_param("tmax") * ns)
{
  // G4 seems to have issues in the um range
  zmin -= copysign(zmin, 1. / 1e6 * cm);
  zmax += copysign(zmax, 1. / 1e6 * cm);
}

PHG4CylinderSteppingAction::~PHG4CylinderSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a nullptr pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4CylinderSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  // get volume of the current step
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4TouchableHandle touchpost = aStep->GetPostStepPoint()->GetTouchableHandle();

  G4VPhysicalVolume* volume = touch->GetVolume();
  // G4 just calls  UserSteppingAction for every step (and we loop then over all our
  // steppingactions. First we have to check if we are actually in our volume
  if (!detector_->IsInCylinder(volume))
  {
    return false;
  }

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;

  const G4Track* aTrack = aStep->GetTrack();

  // if this cylinder stops everything, just put all kinetic energy into edep
  if (IsBlackHole)
  {
    if ((!isfinite(tmin) && !isfinite(tmax)) ||
        aTrack->GetGlobalTime() < tmin ||
        aTrack->GetGlobalTime() > tmax)
    {
      edep = aTrack->GetKineticEnergy() / GeV;
      G4Track* killtrack = const_cast<G4Track*>(aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }
  }

  int layer_id = detector_->get_Layer();
  // test if we are active
  if (active)
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
    //        cout << "time prepoint: " << prePoint->GetGlobalTime()/ns << endl;
    //        cout << "time postpoint: " << postPoint->GetGlobalTime()/ns << endl;
    //        cout << "kinetic energy: " <<  aTrack->GetKineticEnergy()/GeV << endl;
    //       G4ParticleDefinition* def = aTrack->GetDefinition();
    //       cout << "Particle: " << def->GetParticleName() << endl;
    int prepointstatus = prePoint->GetStepStatus();
    if (prepointstatus == fGeomBoundary ||
        prepointstatus == fUndefined ||
        (prepointstatus == fPostStepDoItProc && savepoststepstatus == fGeomBoundary) ||
        use_g4_steps > 0)
    {
      if (prepointstatus == fPostStepDoItProc && savepoststepstatus == fGeomBoundary)
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

      if (!hit)
      {
        hit = new PHG4Hitv1();
      }

      hit->set_layer((unsigned int) layer_id);

      //here we set the entrance values in cm
      hit->set_x(0, prePoint->GetPosition().x() / cm);
      hit->set_y(0, prePoint->GetPosition().y() / cm);
      hit->set_z(0, prePoint->GetPosition().z() / cm);

      hit->set_px(0, prePoint->GetMomentum().x() / GeV);
      hit->set_py(0, prePoint->GetMomentum().y() / GeV);
      hit->set_pz(0, prePoint->GetMomentum().z() / GeV);

      // time in ns
      hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
      //set and save the track ID
      hit->set_trkid(aTrack->GetTrackID());
      savetrackid = aTrack->GetTrackID();
      //set the initial energy deposit
      hit->set_edep(0);
      if (!geantino && !IsBlackHole)
      {
        hit->set_eion(0);
      }
      if (save_light_yield)
      {
        hit->set_light_yield(0);
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

      if (hit->get_z(0) * cm > zmax || hit->get_z(0) * cm < zmin)
      {
	boost::io::ios_precision_saver ips(cout);
        cout << detector_->SuperDetector() << std::setprecision(9)
             << "PHG4CylinderSteppingAction: Entry hit z " << hit->get_z(0) * cm
             << " outside acceptance,  zmin " << zmin
             << ", zmax " << zmax << ", layer: " << layer_id << endl;
      }
    }
    // here we just update the exit values, it will be overwritten
    // for every step until we leave the volume or the particle
    // ceases to exist
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
      cout << "hits do not belong to the same track" << endl;
      cout << "saved track: " << savetrackid
           << ", current trackid: " << aTrack->GetTrackID()
           << endl;
      exit(1);
    }
    saveprestepstatus = prePoint->GetStepStatus();
    savepoststepstatus = postPoint->GetStepStatus();
    savevolpre = volume;
    savevolpost = touchpost->GetVolume();
 
    hit->set_x(1, postPoint->GetPosition().x() / cm);
    hit->set_y(1, postPoint->GetPosition().y() / cm);
    hit->set_z(1, postPoint->GetPosition().z() / cm);

    hit->set_px(1, postPoint->GetMomentum().x() / GeV);
    hit->set_py(1, postPoint->GetMomentum().y() / GeV);
    hit->set_pz(1, postPoint->GetMomentum().z() / GeV);

    hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
    //sum up the energy to get total deposited
    hit->set_edep(hit->get_edep() + edep);
    if (hit->get_z(1) * cm > zmax || hit->get_z(1) * cm < zmin)
    {
      cout << detector_->SuperDetector() << std::setprecision(9)
           << " PHG4CylinderSteppingAction: Exit hit z " << hit->get_z(1) * cm
           << " outside acceptance zmin " << zmin
           << ", zmax " << zmax << ", layer: " << layer_id << endl;
    }
    if (geantino)
    {
      hit->set_edep(-1);  // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
    }
    else
    {
      if (!IsBlackHole)
      {
        double eion = edep - aStep->GetNonIonizingEnergyDeposit() / GeV;
        hit->set_eion(hit->get_eion() + eion);
      }
    }
    if (save_light_yield)
    {
      double light_yield = GetVisibleEnergyDeposition(aStep);
      hit->set_light_yield(hit->get_light_yield() + light_yield);
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
        aTrack->GetTrackStatus() == fStopAndKill ||
        use_g4_steps > 0)
    {
      // save only hits with energy deposit (or -1 for geantino)
      if (hit->get_edep())
      {
        hits_->AddHit(layer_id, hit);
        if (saveshower)
        {
          saveshower->add_g4hit_id(hits_->GetID(), hit->get_hit_id());
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
void PHG4CylinderSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  string hitnodename;
  if (detector_->SuperDetector() != "NONE")
  {
    hitnodename = "G4HIT_" + detector_->SuperDetector();
  }
  else
  {
    hitnodename = "G4HIT_" + detector_->GetName();
  }

  //now look for the map and grab a pointer to it.
  hits_ = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());

  // if we do not find the node we need to make it.
  if (!hits_ && !IsBlackHole)
  {
    std::cout << "PHG4CylinderSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
}
