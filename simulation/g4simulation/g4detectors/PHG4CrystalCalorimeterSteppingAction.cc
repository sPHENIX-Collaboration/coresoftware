#include "PHG4CrystalCalorimeterSteppingAction.h"
#include "PHG4CrystalCalorimeterDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <fun4all/getClass.h>

#include <Geant4/G4Step.hh>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp> // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 4 && BOOST_VERSION == 105700 )
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma message "ignoring bogus gcc warning in boost header lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#pragma GCC diagnostic warning "-Wuninitialized"
#else
#include <boost/lexical_cast.hpp>
#endif

#include <iostream>

using namespace std;


//____________________________________________________________________________..
PHG4CrystalCalorimeterSteppingAction::PHG4CrystalCalorimeterSteppingAction( PHG4CrystalCalorimeterDetector* detector ):
  PHG4SteppingAction(NULL),
  detector_( detector ),
  hits_(NULL),
  absorberhits_(NULL),
  hit(NULL)
{

}


//____________________________________________________________________________..
bool PHG4CrystalCalorimeterSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch->GetVolume();

  // detector_->IsInCrystalCalorimeter(volume)
  // returns
  //  0 is outside of Crystal Calorimeter
  //  1 is inside scintillator scrystal
  // -1 is absorber (dead material)

  int whichactive = detector_->IsInCrystalCalorimeter(volume);

  if ( !whichactive  )
    {
      return false;
    }

  int layer_id = detector_->get_Layer();
  int tower_id = -1;
  int idx_j = -1;
  int idx_k = -1;
  int idx_l = -1;

  if (whichactive > 0) // in crystal
    {
      /* Find indizes of crystal containing this step */
      /** @TODO ParseName works for planar geometry, WhatAreYou works for projective geometry
       * (because of the additional layer of logical volume). Need to clean this up. */
      //WhatAreYou(touch, idx_j, idx_k);
      ParseName(touch->GetVolume(1), idx_j, idx_k);
      tower_id = touch->GetCopyNumber();
    }
  else if (whichactive < 0)
    {
      //no tower ID for absorber
      tower_id = touch->GetCopyNumber();
    }
  else
    {
      tower_id = touch->GetCopyNumber();
    }

  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  /* Get pointer to associated Geant4 track */
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (detector_->IsBlackHole())
    {
      edep = aTrack->GetKineticEnergy() / GeV;
      G4Track* killtrack = const_cast<G4Track *> (aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }

  /* Make sure we are in a volume */
  if ( detector_->IsActive() )
    {
      /* Check if particle is 'geantino' */
      bool geantino = false;
      if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
	  aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
	{
	  geantino = true;
	}

      /* Get Geant4 pre- and post-step points */
      G4StepPoint * prePoint = aStep->GetPreStepPoint();
      G4StepPoint * postPoint = aStep->GetPostStepPoint();

      switch (prePoint->GetStepStatus())
	{
	case fGeomBoundary:
	case fUndefined:
	  hit = new PHG4Hitv1();
//	  hit->set_layer(0);
	  hit->set_scint_id(tower_id);

	  /* Set hit location (tower index) only for crystals*/
	  if (whichactive > 0)
	    {
	      hit->set_index_j(idx_j);
	      hit->set_index_k(idx_k);
	      hit->set_index_l(idx_l);
	    }

	  /* Set hit location (space point) only for crystals*/
	  hit->set_x( 0, prePoint->GetPosition().x() / cm);
	  hit->set_y( 0, prePoint->GetPosition().y() / cm );
	  hit->set_z( 0, prePoint->GetPosition().z() / cm );

	  /* Set hit time */
	  hit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );

	  /* set the track ID */
	  {
	    int trkoffset = 0;
	    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	      {
		if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		  {
		    trkoffset = pp->GetTrackIdOffset();
		  }
	      }
	    hit->set_trkid(aTrack->GetTrackID() + trkoffset);
	  }

	  /* set intial energy deposit */
	  hit->set_edep( 0 );
	  hit->set_eion( 0 );

	  /* Now add the hit to the hit collection */
	  if (whichactive == 1) // return of IsInInnerHcalDetector, > 0 hit in scintillator, < 0 hit in absorber
	    {
	      // Now add the hit
	      hits_->AddHit(layer_id, hit);
	    }
	  else
	    {
	      absorberhits_->AddHit(layer_id, hit);
	    }
	  break;
	default:
	  break;
	}

      /* Update exit values- will be overwritten with every step until
       * we leave the volume or the particle ceases to exist */
      hit->set_x( 1, postPoint->GetPosition().x() / cm );
      hit->set_y( 1, postPoint->GetPosition().y() / cm );
      hit->set_z( 1, postPoint->GetPosition().z() / cm );

      hit->set_t( 1, postPoint->GetGlobalTime() / nanosecond );

      /* sum up the energy to get total deposited */
      hit->set_edep(hit->get_edep() + edep);
      hit->set_eion(hit->get_eion() + eion);

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

      return true;

    }
  else
    {
      return false;
    }
}


//____________________________________________________________________________..
void PHG4CrystalCalorimeterSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
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

  // if we do not find the node it's messed up.
  if ( ! hits_ )
    {
      std::cout << "PHG4CrystalCalorimeterSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
    }
  if ( ! absorberhits_)
    {
      if (verbosity > 0)
	{
	  cout << "PHG4CrystalCalorimeterSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
	}
    }

}

int
PHG4CrystalCalorimeterSteppingAction::WhatAreYou(G4TouchableHandle touch, int& j, int& k)
{
        int j_0, k_0;           //The k and k indices for the crystal within the 2x2 matrix
        int j_1, k_1;           //The k and k indices for the 2x2 within the 4x4 
        int j_2, k_2;           //The k and k indices for the 4x4 within the mother volume
        //int j, k;             //The final indices of the crystal

	string name = touch->GetVolume(0)->GetName();
	string name_2 = touch->GetVolume(1)->GetName();

	if (name.find("rystal") != string::npos)
	{
		G4VPhysicalVolume* crystal = touch->GetVolume(0);		//Get the crystal solid
		G4VPhysicalVolume* TwoByTwo = touch->GetVolume(1);		//Get the crystal solid
		G4VPhysicalVolume* FourByFour = touch->GetVolume(2);		//Get the crystal solid

		ParseName(crystal, j_0, k_0);
		ParseName(TwoByTwo, j_1, k_1);
		ParseName(FourByFour, j_2, k_2);
		
		j = (j_0*1) + (j_1*2) + (j_2*4);
		k = (k_0*1) + (k_1*2) + (k_2*4);

	}
	else if ( (name.find("arbon") != string::npos) )
	{
		G4VPhysicalVolume* carbon = touch->GetVolume(1);		//Get the carbon fiber solid
		ParseName(carbon, j_0, k_0);
		j = j_0;
		k = k_0;
		//cout << "Carbon_" << j << "_" << k << endl;
	}
	else
	{
		cout << "ERROR!!! RUN FOR YOUR LIVES!!!" << endl;
		j = -1;
		k = -1;
	}

        return 0;
}

int 
PHG4CrystalCalorimeterSteppingAction::ParseName( G4VPhysicalVolume* volume, int& j, int& k ) 
{
	boost::char_separator<char> sep("_");
	boost::tokenizer<boost::char_separator<char> > tok(volume->GetName(), sep);
	boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
	for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
	{
		if (*tokeniter == "j")
		{
			++tokeniter;
			j = boost::lexical_cast<int>(*tokeniter);
		}
		else if (*tokeniter == "k")
		{
			++tokeniter;
			k = boost::lexical_cast<int>(*tokeniter);
		}
	}
	
	return 0;
}
