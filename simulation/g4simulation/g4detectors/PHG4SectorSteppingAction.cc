#include "PHG4SectorSteppingAction.h"
#include "PHG4SectorDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <fun4all/getClass.h>

#include <Geant4/G4Step.hh>

#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4SectorSteppingAction::PHG4SectorSteppingAction(PHG4SectorDetector* detector) :
    detector_(detector)
{
}

//____________________________________________________________________________..
bool
PHG4SectorSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{

  // get volume of the current step
  G4VPhysicalVolume* volume =
      aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit()-aStep->GetNonIonizingEnergyDeposit()) / GeV;

  const G4Track* aTrack = aStep->GetTrack();

  int layer_id = 0;
  // make sure we are in a volume
  if (detector_->IsInSectorActive(volume))
    {
      bool geantino = false;
      // the check for the pdg code speeds things up, I do not want to make 
      // an expensive string compare for every track when we know
      // geantino or chargedgeantino has pid=0
      if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0
          && aTrack->GetParticleDefinition()->GetParticleName().find("geantino")
              != string::npos)
        {
          geantino = true;
        }
      G4StepPoint * prePoint = aStep->GetPreStepPoint();
      G4StepPoint * postPoint = aStep->GetPostStepPoint();
//       cout << "track id " << aTrack->GetTrackID() << endl;
//       cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
//       cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
      switch (prePoint->GetStepStatus())
        {
      case fGeomBoundary:
      case fUndefined:
        hit = new PHG4Hitv1();
        //here we set the entrance values in cm
        hit->set_x(0, prePoint->GetPosition().x() / cm);
        hit->set_y(0, prePoint->GetPosition().y() / cm);
        hit->set_z(0, prePoint->GetPosition().z() / cm);
        // time in ns
        hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
        //set the track ID
          {
            int trkoffset = 0;
            if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
              {
                if (PHG4TrackUserInfoV1* pp =
                    dynamic_cast<PHG4TrackUserInfoV1*>(p))
                  {
                    trkoffset = pp->GetTrackIdOffset();
                  }
              }
            hit->set_trkid(aTrack->GetTrackID() + trkoffset);
          }

        //set the initial energy deposit
        hit->set_edep(0);
        hit->set_eion(0); // only implemented for v5 otherwise empty
//        hit->set_light_yield(0);

        //layer_id is sector number
        layer_id = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1);

        // Now add the hit
        hits_->AddHit(layer_id, hit);

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
      //sum up the energy to get total deposited
      hit->set_edep(hit->get_edep() + edep);
      hit->set_eion(hit->get_eion() + eion);
      if (geantino)
        {
          hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
        }
      if (edep > 0)
        {
          if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
            {
              if (PHG4TrackUserInfoV1* pp =
                  dynamic_cast<PHG4TrackUserInfoV1*>(p))
                {
                  pp->SetKeep(1); // we want to keep the track
                }
            }
        }
      hit->set_path_length(aTrack->GetTrackLength()/cm);

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
void
PHG4SectorSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  if (!hits_)
    {
      std::cout << "PHG4SectorSteppingAction::SetTopNode - unable to find "
          << hitnodename << std::endl;
    }

}
