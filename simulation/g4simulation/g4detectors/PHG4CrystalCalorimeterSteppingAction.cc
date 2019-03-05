#include "PHG4CrystalCalorimeterSteppingAction.h"
#include "PHG4CrystalCalorimeterDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4MaterialCutsCouple.hh>
#include <Geant4/G4SystemOfUnits.hh>

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
  detector_( detector ),
  hits_(NULL),
  absorberhits_(NULL),
  hit(NULL),
  savehitcontainer(NULL),
  saveshower(NULL)
{

}

PHG4CrystalCalorimeterSteppingAction::~PHG4CrystalCalorimeterSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a NULL pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
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

  if (whichactive > 0) // in crystal
    {
      /* Find indizes of crystal containing this step */
      if ( touch->GetVolume(2)->GetName().find("_j_") != string::npos )
	FindTowerIndex2LevelUp(touch, idx_j, idx_k);
      else
	FindTowerIndex(touch, idx_j, idx_k);

      tower_id = touch->GetCopyNumber();
    }
  else
    {
      tower_id = touch->GetCopyNumber();
    }

  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;
  G4double light_yield = 0;

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
      int idx_l = -1;
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
	  if (! hit)
	    {
	      hit = new PHG4Hitv1();
	    }
	  hit->set_scint_id(tower_id);

	  /* Set hit location (tower index) */
	  hit->set_index_j(idx_j);
	  hit->set_index_k(idx_k);
	  hit->set_index_l(idx_l);

	  /* Set hit location (space point)*/
	  hit->set_x( 0, prePoint->GetPosition().x() / cm);
	  hit->set_y( 0, prePoint->GetPosition().y() / cm );
	  hit->set_z( 0, prePoint->GetPosition().z() / cm );

	  /* Set hit time */
	  hit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );

	  /* Set the track ID */
	  hit->set_trkid(aTrack->GetTrackID());

	  /* Set intial energy deposit */
	  hit->set_edep( 0 );
	  hit->set_eion( 0 );

	  /* Now add the hit to the hit collection */
	  // here we do things which are different between scintillator and absorber hits
	  if (whichactive > 0)
	    {
              savehitcontainer = hits_;
	      hit->set_light_yield( 0 ); // for scintillator only, initialize light yields
	    }
	  else
	    {
	      savehitcontainer = absorberhits_;
	    }

	  // here we set what is common for scintillator and absorber hits
	  if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	    {
	      if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		{
		  hit->set_trkid(pp->GetUserTrackId());
		  hit->set_shower_id(pp->GetShower()->get_id());
		  saveshower = pp->GetShower();
		}
	    }
	  break;
	default:
	  break;
	}

      if (whichactive > 0)
	{
	  /* Use 'light yield = depoisted energy (ionization)' for this calorimeter.
	   * At this point, the calorimeters using orgqanic scintillators apply Birk's law correction via
	   * light_yield = GetVisibleEnergyDeposition(aStep);
	   * to account for its effect.
	   */
	  light_yield = eion;
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
      if (whichactive > 0)
	{
	  hit->set_light_yield(hit->get_light_yield() + light_yield);
	}

      if (geantino)
	{
	  hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
          hit->set_eion(-1);
	  if (whichactive > 0)
	    {
	      hit->set_light_yield(-1);
	    }
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
      // if any of these conditions is true this is the last step in
      // this volume and we need to save the hit
      // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
      // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
      // (not sure if this will ever be the case)
      // aTrack->GetTrackStatus() == fStopAndKill: track ends
      if (postPoint->GetStepStatus() == fGeomBoundary || postPoint->GetStepStatus() == fWorldBoundary|| aTrack->GetTrackStatus() == fStopAndKill)
	{
          // save only hits with energy deposit (or -1 for geantino)
	  if (hit->get_edep())
	    {
	      savehitcontainer->AddHit(layer_id, hit);
	      if (saveshower)
		{
		  saveshower->add_g4hit_id(hits_->GetID(),hit->get_hit_id());
		}
	      // ownership has been transferred to container, set to null
	      // so we will create a new hit for the next track
	      hit = NULL;
	    }
	  else
	    {
	      // if this hit has no energy deposit, just reset it for reuse
	      // this means we have to delete it in the dtor. If this was
	      // the last hit we processed the memory is still allocated
	      hit->Reset();
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
      if (Verbosity() > 0)
	{
	  cout << "PHG4CrystalCalorimeterSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
	}
    }

}

int
PHG4CrystalCalorimeterSteppingAction::FindTowerIndex(G4TouchableHandle touch, int& j, int& k)
{
        int j_0, k_0;           //The k and k indices for the scintillator / tower

        G4VPhysicalVolume* tower = touch->GetVolume(1);		//Get the tower solid
	ParseG4VolumeName(tower, j_0, k_0);

        j = (j_0*1);
        k = (k_0*1);

        return 0;
}

int
PHG4CrystalCalorimeterSteppingAction::FindTowerIndex2LevelUp(G4TouchableHandle touch, int& j, int& k)
{
        int j_0, k_0;           //The k and k indices for the crystal within the 2x2 matrix

	string name = touch->GetVolume(0)->GetName();

	if (name.find("rystal") != string::npos)
	{
		G4VPhysicalVolume* crystal = touch->GetVolume(0);		//Get the crystal solid
		G4VPhysicalVolume* TwoByTwo = touch->GetVolume(1);		//Get the crystal solid
		G4VPhysicalVolume* FourByFour = touch->GetVolume(2);		//Get the crystal solid
                int j_1, k_1;  //The k and k indices for the 2x2 within the 4x4
                int j_2, k_2;  //The k and k indices for the 4x4 within the mother volume
		ParseG4VolumeName(crystal, j_0, k_0);
		ParseG4VolumeName(TwoByTwo, j_1, k_1);
		ParseG4VolumeName(FourByFour, j_2, k_2);

		j = (j_0*1) + (j_1*2) + (j_2*4);
		k = (k_0*1) + (k_1*2) + (k_2*4);

	}
	else if ( (name.find("arbon") != string::npos) )
	{
		G4VPhysicalVolume* carbon = touch->GetVolume(1);		//Get the carbon fiber solid
		ParseG4VolumeName(carbon, j_0, k_0);
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
PHG4CrystalCalorimeterSteppingAction::ParseG4VolumeName( G4VPhysicalVolume* volume, int& j, int& k )
{
	boost::char_separator<char> sep("_");
	boost::tokenizer<boost::char_separator<char> > tok(volume->GetName(), sep);
	boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
	for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
	{
		if (*tokeniter == "j")
		{
			++tokeniter;
			if (tokeniter == tok.end()) break;
			j = boost::lexical_cast<int>(*tokeniter);
		}
		else if (*tokeniter == "k")
		{
			++tokeniter;
      if (tokeniter == tok.end()) break;
			k = boost::lexical_cast<int>(*tokeniter);
		}
	}

	return 0;
}
