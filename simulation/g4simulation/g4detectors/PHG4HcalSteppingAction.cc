#include "PHG4HcalSteppingAction.h"
#include "PHG4HcalDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4HitDefs.h>

#include <g4main/PHG4TrackUserInfoV1.h>


#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4MaterialCutsCouple.hh>

#include <iostream>

using namespace std;
//____________________________________________________________________________..
PHG4HcalSteppingAction::PHG4HcalSteppingAction( PHG4HcalDetector* detector ):
  PHG4SteppingAction(0),
  detector_( detector ),
  hits_(nullptr),
  absorberhits_(nullptr),
  hit(nullptr),
  zmin(NAN),
  zmax(NAN),
  light_scint_model_(true),
  light_balance_(false),
  light_balance_inner_radius_(0.0),
  light_balance_inner_corr_(1.0),
  light_balance_outer_radius_(10.0),
  light_balance_outer_corr_(1.0)  
{}

//____________________________________________________________________________..
bool PHG4HcalSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{

  // get volume of the current step
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit()-aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;

  const G4Track* aTrack = aStep->GetTrack();

  int layer_id = detector_->get_Layer();
  // make sure we are in a volume
  // IsInCylinderActive returns the number of the scintillator 
  // slat which has fired
  int isactive = detector_->IsInCylinderActive(volume);
  if ( isactive > PHG4HcalDetector::INACTIVE ) 
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
      G4StepPoint * prePoint = aStep->GetPreStepPoint();
      G4StepPoint * postPoint = aStep->GetPostStepPoint();
      //       cout << "track id " << aTrack->GetTrackID() << endl;
      //        cout << "time prepoint: " << prePoint->GetGlobalTime() << endl;
      //        cout << "time postpoint: " << postPoint->GetGlobalTime() << endl;
      switch (prePoint->GetStepStatus())
        {
        case fGeomBoundary:
        case fUndefined:

	  hit = new PHG4Hitv1();

	  hit->set_layer((unsigned int)layer_id);
	  hit->set_scint_id(isactive); // isactive contains the scintillator slat id
          //here we set the entrance values in cm
          hit->set_x( 0, prePoint->GetPosition().x() / cm );
          hit->set_y( 0, prePoint->GetPosition().y() / cm );
          hit->set_z( 0, prePoint->GetPosition().z() / cm );

	  // time in ns
          hit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );
	  //set the track ID
	  {
            hit->set_trkid(aTrack->GetTrackID());
            if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	      {
		if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		  {
		    hit->set_trkid(pp->GetUserTrackId());
		    hit->set_shower_id(pp->GetShower()->get_id());
		  }
	      }
	  }
          //set the initial energy deposit
          hit->set_edep(0);
	  hit->set_eion(0); // only implemented for v5 otherwise empty
	  hit->set_light_yield(0);
	  //	  hit->print();
          // Now add the hit
	  if (isactive >= 0) // the slat ids start with zero
	    {
	      //	      unsigned int shift_layer_id = layer_id << (phg4hitdefs::keybits - 3);
	      hits_->AddHit(layer_id, hit); // scintillator id is coded into layer number

	      {
		if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
		  {
		    if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		      {
			pp->GetShower()->add_g4hit_id(hits_->GetID(),hit->get_hit_id());
		      }
		  }
	      }
	      
	    }
	  else
	    {
	      absorberhits_->AddHit(layer_id, hit);

	      {
		if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
		  {
		    if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		      {
			pp->GetShower()->add_g4hit_id(absorberhits_->GetID(),hit->get_hit_id());
		      }
		  }
	      }
	    }

	  if (hit->get_z(0) > zmax || hit->get_z(0) < zmin)
	    {
	      cout << "PHG4HcalSteppingAction: hit outside acceptance, layer: " << layer_id << endl;
	      hit->identify();
	    }
          break;
        default:
          break;
        }
      // here we just update the exit values, it will be overwritten
      // for every step until we leave the volume or the particle
      // ceases to exist
      hit->set_x( 1, postPoint->GetPosition().x() / cm );
      hit->set_y( 1, postPoint->GetPosition().y() / cm );
      hit->set_z( 1, postPoint->GetPosition().z() / cm );

      hit->set_t( 1, postPoint->GetGlobalTime() / nanosecond );

      if (isactive >= 0) // the slat ids start with zero
        {

          if (light_scint_model_)
            {
              light_yield = GetVisibleEnergyDeposition(aStep);

              static bool once = true;
              if (once && edep>0)
                {
                  once = false;

                  cout << "PHG4HcalSteppingAction::UserSteppingAction::"
		    //
		       << detector_->GetName() << " - "
		       << " use scintillating light model at each Geant4 steps. "
		       <<"First step: "
		       <<"Material = "<<aTrack->GetMaterialCutsCouple()->GetMaterial()->GetName()<<", "
		       <<"Birk Constant = "<<aTrack->GetMaterialCutsCouple()->GetMaterial()->GetIonisation()->GetBirksConstant()<<","
		       <<"edep = " <<edep<<", "
		       <<"eion = " <<eion<<", "
		       <<"light_yield = " <<light_yield
		       << endl;
                }

            }
          else
            {
              light_yield = eion;
            }

          if (light_balance_)
            {
              float r = sqrt(
			     pow(postPoint->GetPosition().x() / cm, 2)
			     + pow(postPoint->GetPosition().y() / cm, 2));
              const float cor = GetLightCorrection(r);
              light_yield = light_yield * cor;

              static bool once = true;
              if (once && light_yield>0)
                {
                  once = false;

                  cout << "PHG4HcalSteppingAction::UserSteppingAction::"
		    //
		       << detector_->GetName() << " - "
		       << " use a simple light collection model with linear radial dependence. "
		       <<"First step: "
		       <<"r = " <<r<<", "
		       <<"correction ratio = " <<cor<<", "
		       <<"light_yield after cor. = " <<light_yield
		       << endl;
                }

            }
        }
      
      //sum up the energy to get total deposited
      hit->set_edep(hit->get_edep() + edep);
      hit->set_eion(hit->get_eion() + eion);
      hit->set_light_yield(hit->get_light_yield() + light_yield);
      hit->set_path_length(aTrack->GetTrackLength() / cm);

      if (hit->get_z(1) > zmax || hit->get_z(1) < zmin)
	{
	  cout << "PHG4HcalSteppingAction: hit outside acceptance zmin " << zmin << ", zmax " << zmax << " at exit" << endl;
	  hit->identify();
	}
      if (geantino)
	{
	  hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
          hit->set_eion(-1);
	}
      if (edep > 0)
	{
	  if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	    {
	      if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		{
		  pp->SetKeep(1); // we want to keep the track
		}
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
void PHG4HcalSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
{

  string hitnodename;
  string absorbernodename;
  if (detector_->SuperDetector() != "NONE")
    {
      hitnodename = "G4HIT_" + detector_->SuperDetector();
      absorbernodename =  "G4HIT_ABSORBER_" + detector_->SuperDetector();
    }
  else
    {
      hitnodename = "G4HIT_" + detector_->GetName();
      absorbernodename =  "G4HIT_ABSORBER_" + detector_->GetName();
    }

  //now look for the map and grab a pointer to it.
  hits_ =  findNode::getClass<PHG4HitContainer>( topNode , hitnodename.c_str() );
  absorberhits_ =  findNode::getClass<PHG4HitContainer>( topNode , absorbernodename.c_str() );
  // if we do not find the node we need to make it.
  if ( ! hits_ )
    {
      std::cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
    }
  if ( ! absorberhits_)
    {
      if (Verbosity() > 0)
	{
	  std::cout << "PHG4HcalSteppingAction::SetTopNode - unable to find " << absorbernodename << std::endl;
	}
    }
}

float PHG4HcalSteppingAction::GetLightCorrection(float r) {
  float m = (light_balance_outer_corr_ - light_balance_inner_corr_)/(light_balance_outer_radius_ - light_balance_inner_radius_);
  float b = light_balance_inner_corr_ - m*light_balance_inner_radius_;
  float value = m*r+b;  
  if (value > 1.0) return 1.0;
  if (value < 0.0) return 0.0;

  return value;
}
