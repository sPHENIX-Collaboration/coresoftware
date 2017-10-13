/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.6 $$
 * \date $$Date: 2015/01/07 23:50:05 $$
 */

#include "PHG4SpacalPrototypeSteppingAction.h"
#include "PHG4SpacalPrototypeDetector.h"
#include "PHG4CylinderGeom_Spacalv3.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4HitDefs.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4SystemOfUnits.hh>

#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4SpacalPrototypeSteppingAction::PHG4SpacalPrototypeSteppingAction(PHG4SpacalPrototypeDetector* detector) :
  PHG4SteppingAction(0), 
  detector_(detector), 
  hits_(nullptr), 
  absorberhits_(nullptr),
  hit(nullptr),
  savehitcontainer(nullptr),
  saveshower(nullptr)
{}

//____________________________________________________________________________..
bool
PHG4SpacalPrototypeSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{

  // get volume of the current step
  G4VPhysicalVolume* volume =
      aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit()
      - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  const G4Track* aTrack = aStep->GetTrack();

  const int layer_id = 0;
  // make sure we are in a volume
  // IsInCylinderActive returns the number of the scintillator 
  // slat which has fired
  int isactive = detector_->IsInCylinderActive(volume);
  if (isactive > PHG4SpacalPrototypeDetector::INACTIVE)
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
      int scint_id = -1;

      if (//
          detector_->get_geom()->get_config() == PHG4SpacalPrototypeDetector::SpacalGeom_t::kFullProjective_2DTaper //
          or //
      detector_->get_geom()->get_config() == PHG4SpacalPrototypeDetector::SpacalGeom_t::kFullProjective_2DTaper_SameLengthFiberPerTower//
      )
        {
          //SPACAL ID that is associated with towers
          int sector_ID =0;
          int tower_ID = 0;
          int fiber_ID = 0;

          if (isactive == PHG4SpacalPrototypeDetector::FIBER_CORE)
            {

              fiber_ID = prePoint->GetTouchable()->GetReplicaNumber(1);
              tower_ID = prePoint->GetTouchable()->GetReplicaNumber(2);
              sector_ID  = prePoint->GetTouchable()->GetReplicaNumber(3);

            }

          else if (isactive == PHG4SpacalPrototypeDetector::FIBER_CLADING)
            {
              fiber_ID = prePoint->GetTouchable()->GetReplicaNumber(0);
              tower_ID = prePoint->GetTouchable()->GetReplicaNumber(1);
              sector_ID  = prePoint->GetTouchable()->GetReplicaNumber(2);
            }

          else if (isactive == PHG4SpacalPrototypeDetector::ABSORBER)
            {
              tower_ID = prePoint->GetTouchable()->GetReplicaNumber(0);
              sector_ID  = prePoint->GetTouchable()->GetReplicaNumber(1);
            }

          // compact the tower/sector/fiber ID into 32 bit scint_id, so we could save some space for SPACAL hits
          scint_id = PHG4CylinderGeom_Spacalv3::scint_id_coder(sector_ID, tower_ID, fiber_ID).scint_ID;

        }
      else
        {
          // other configuraitons
          if (isactive == PHG4SpacalPrototypeDetector::FIBER_CORE)
	    {
	      scint_id = prePoint->GetTouchable()->GetReplicaNumber(2);
	    }
          else if (isactive == PHG4SpacalPrototypeDetector::FIBER_CLADING)
	    {
	      scint_id = prePoint->GetTouchable()->GetReplicaNumber(1);
	    }
          else
	    {
	      scint_id = prePoint->GetTouchable()->GetReplicaNumber(0);
	    }
        }

      //       cout << "track id " << aTrack->GetTrackID() << endl;
      //        cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
      //        cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
      switch (prePoint->GetStepStatus())
        {
        case fGeomBoundary:
        case fUndefined:
	  //	case fPostStepDoItProc: // from point like interaction (compton,decay,hard rad)

	if (! hit)
	  {
	    hit = new PHG4Hitv1();
	  }
	hit->set_layer((unsigned int) layer_id);
        hit->set_scint_id(scint_id); // isactive contains the scintillator slat id
        //here we set the entrance values in cm
        hit->set_x(0, prePoint->GetPosition().x() / cm);
        hit->set_y(0, prePoint->GetPosition().y() / cm);
        hit->set_z(0, prePoint->GetPosition().z() / cm);

        // time in ns
        hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);

        if (isactive == PHG4SpacalPrototypeDetector::FIBER_CORE) // only for active areas
          {
            // store all pre local coordinates
            StoreLocalCoordinate(hit, aStep, true, false);
          }

	//set the track ID
	  hit->set_trkid(aTrack->GetTrackID());
	  if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	    {
	      if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		{
		  hit->set_trkid(pp->GetUserTrackId());
		  hit->set_shower_id(pp->GetShower()->get_id());
		}
	    }
        //set the initial energy deposit
        hit->set_edep(0);
        if (isactive == PHG4SpacalPrototypeDetector::FIBER_CORE) // only for active areas
          {
            hit->set_eion(0); // only for fiber core hits
            hit->set_light_yield(0);
	    savehitcontainer = hits_;
          }
	else
	  {
	    savehitcontainer = absorberhits_;
          }

	if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	  {
	    if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
	      {
		hit->set_trkid(pp->GetUserTrackId());
		hit->set_shower_id(pp->GetShower()->get_id());
		saveshower =  pp->GetShower();
	      }
	  }

//        if (hit->get_z(0) > get_zmax() || hit->get_z(0) < get_zmin())
//          {
//            cout << "PHG4SpacalPrototypeSteppingAction: hit outside acceptance, layer: "
//                << layer_id << endl;
//            hit->identify();
//          }
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

      if (isactive == PHG4SpacalPrototypeDetector::FIBER_CORE) // only for active areas
        {
          // store all pre local coordinates
          StoreLocalCoordinate(hit, aStep, false, true);
          hit->set_eion(hit->get_eion() + eion);

          double light_yield = GetVisibleEnergyDeposition(aStep);

          static bool once = true;
          if (once and edep > 0)
            {
              once = false;

	      if (verbosity > 0) {
		cout << "PHG4SpacalPrototypeSteppingAction::UserSteppingAction::"
                  //
		     << detector_->GetName() << " - "
		     << " use scintillating light model at each Geant4 steps. "
		     << "First step: " << "Material = "
		     << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetName()
		     << ", " << "Birk Constant = "
		     << aTrack->GetMaterialCutsCouple()->GetMaterial()->GetIonisation()->GetBirksConstant()
		     << "," << "edep = " << edep << ", " << "eion = " << eion
		     << ", " << "light_yield = " << light_yield << endl;
	      }
	    }

          hit->set_light_yield(hit->get_light_yield() + light_yield);
        }

//      if (hit->get_z(1) > get_zmax() || hit->get_z(1) < get_zmin())
//        {
//          cout << "PHG4SpacalPrototypeSteppingAction: hit outside acceptance get_zmin() "
//              << get_zmin() << ", get_zmax() " << get_zmax() << " at exit"
//              << endl;
//          hit->identify();
//        }
      if (geantino)
        {
          hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
//          hit->set_eion(-1);
        }
      if (edep > 0)
        {
          if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
            {
              if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
                {
                  pp->SetKeep(1); // we want to keep the track
                }
            }
        }

      // if any of these conditions is true this is the last step in
      // this volume and we need to save the hit
      // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
      // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
      // (not sure if this will ever be the case)
      // aTrack->GetTrackStatus() == fStopAndKill: track ends
      if (postPoint->GetStepStatus() == fGeomBoundary || 
          postPoint->GetStepStatus() == fWorldBoundary|| 
          aTrack->GetTrackStatus() == fStopAndKill)
	{
          // save only hits with energy deposit (or -1 for geantino)
	  if (hit->get_edep())
	    {
	      savehitcontainer->AddHit(layer_id, hit);
	      if (saveshower)
		{
		  // save shower under container id of active volume (for running without absorber hit container to save space)
		  saveshower->add_g4hit_id(savehitcontainer->GetID(),hit->get_hit_id());
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
void
PHG4SpacalPrototypeSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
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
  absorberhits_ = findNode::getClass<PHG4HitContainer>(topNode,
      absorbernodename.c_str());
  // if we do not find the node we need to make it.
  if (!hits_)
    {
      std::cout << "PHG4SpacalPrototypeSteppingAction::SetTopNode - unable to find "
          << hitnodename << std::endl;
    }
  if (!absorberhits_)
    {
      if (verbosity > 1)
        {
          std::cout << "PHG4SpacalPrototypeSteppingAction::SetTopNode - unable to find "
              << absorbernodename << std::endl;
        }
    }
}

double
PHG4SpacalPrototypeSteppingAction::get_zmin()
{
  if (!detector_)
    return 0;
  else
    return detector_->get_geom()->get_zmin() - .0001;
}

double
PHG4SpacalPrototypeSteppingAction::get_zmax()
{
  if (!detector_)
    return 0;
  else
    return detector_->get_geom()->get_zmax() + .0001;
}
