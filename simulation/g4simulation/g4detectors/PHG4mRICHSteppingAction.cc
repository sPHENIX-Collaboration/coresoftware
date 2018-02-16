/*===============================================================*
 *                       March 19th 2017                         *
    mRICH Stepping Action created by Cheuk-Ping Wong @GSU        *
 *===============================================================*
 *                       March 19th 2017                         * 
 *---------------------------------------------------------------*
 * Modified from PHG4ForwardEcalSteppingAction.cc This code still*
 * need polishing                                                *
 *===============================================================*/
#include "PHG4mRICHSteppingAction.h"
#include "PHG4mRICHDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4MaterialCutsCouple.hh>

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
using namespace CLHEP;

//____________________________________________________________________________..
PHG4mRICHSteppingAction::PHG4mRICHSteppingAction( PHG4mRICHDetector* detector,PHParameters* params):
  detector_( detector ),
  active(params->get_int_param("active")),
  IsBlackHole(params->get_int_param("blackhole")),
  use_g4_steps(params->get_int_param("use_g4steps")),
  detectorname(params->get_string_param("detectorname")),
  superdetector(params->get_string_param("superdetector")),
  hits_(NULL),
  absorberhits_(NULL),
  hit(NULL)
{}

//____________________________________________________________________________..
bool PHG4mRICHSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch->GetVolume();

  // detector_->IsInForwardEcal(volume)
  // returns
  //  0 is outside of Forward ECAL
  //  1 is inside scintillator
  // -1 is inside absorber (dead material)

  int whichactive = detector_->IsInmRICH(volume);

  if ( !whichactive  ) return false;

  /* Get energy deposited by this step */
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

  /* Get pointer to associated Geant4 track */
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (IsBlackHole) {
    edep = aTrack->GetKineticEnergy() / GeV;
    G4Track* killtrack = const_cast<G4Track *> (aTrack);
    killtrack->SetTrackStatus(fStopAndKill);
  }

  /* Make sure we are in a volume */
  if ( active ) {
    /* Check if particle is 'geantino' */
    bool geantino = false;
    if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
	aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos) {
      geantino = true;
    }

    /* Get Geant4 pre- and post-step points */
    G4StepPoint * prePoint = aStep->GetPreStepPoint();
    G4StepPoint * postPoint = aStep->GetPostStepPoint();

    switch (prePoint->GetStepStatus())
      {
      //-----------------
      case fGeomBoundary:
      //-----------------
      case fUndefined:
	if (! hit) hit = new PHG4Hitv1();

	/* Set hit location (space point) */
	hit->set_x( 0, prePoint->GetPosition().x() / cm);
	hit->set_y( 0, prePoint->GetPosition().y() / cm );
	hit->set_z( 0, prePoint->GetPosition().z() / cm );
	
	/* Set hit time */
	hit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );
	
	//set the track ID
	hit->set_trkid(aTrack->GetTrackID());
	/* set intial energy deposit */
	hit->set_edep( 0 );
	hit->set_eion( 0 );
	
	// here we set what is common for scintillator and absorber hits
	if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() ) {
	  if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) ) {
	    hit->set_trkid(pp->GetUserTrackId());
	    hit->set_shower_id(pp->GetShower()->get_id());
	  }
	}
	break;
      //-----------------
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
    
    if (geantino) {
      hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
      hit->set_eion(-1);
    }
    if (edep > 0 /*&& (whichactive > 0 || absorbertruth > 0)*/) {
      if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() ) {
	if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) ) {
	  pp->SetKeep(1); // we want to keep the track
	}
      }
    }
    return true;
    
  }
  else return false;

}


//____________________________________________________________________________..
void PHG4mRICHSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
{

  string hitnodename;
  string absorbernodename;

  if (superdetector !="NONE") {
    hitnodename = "G4HIT_" + superdetector;
    absorbernodename =  "G4HIT_ABSORBER_" + superdetector;
  }
  else {
    hitnodename = "G4HIT_" + detectorname;
    absorbernodename =  "G4HIT_ABSORBER_" + detectorname;
  }
  
  //now look for the map and grab a pointer to it.
  hits_ =  findNode::getClass<PHG4HitContainer>( topNode , hitnodename.c_str() );
  absorberhits_ =  findNode::getClass<PHG4HitContainer>( topNode , absorbernodename.c_str() );
  
  // if we do not find the node it's messed up.
  if ( ! hits_ ) {
    std::cout << "PHG4mRICHSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
  }
  if ( ! absorberhits_) {
    if (verbosity > 0) {
      cout << "PHG4mRICHSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
    }
  }
}
//____________________________________________________________________________..
/*
int
PHG4ForwardEcalSteppingAction::FindTowerIndex(G4TouchableHandle touch, int& j, int& k)
{
        int j_0, k_0;           //The j and k indices for the scintillator / tower

        G4VPhysicalVolume* tower = touch->GetVolume(0);		//Get the tower solid
	ParseG4VolumeName(tower, j_0, k_0);

        j = (j_0*1);
        k = (k_0*1);

        return 0;
}

int
PHG4ForwardEcalSteppingAction::ParseG4VolumeName( G4VPhysicalVolume* volume, int& j, int& k ) 
{
	boost::char_separator<char> sep("_");
	boost::tokenizer<boost::char_separator<char> > tok(volume->GetName(), sep);
	boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
	for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
	{
		if (*tokeniter == "j")
		{
			++tokeniter;
			if ( tokeniter == tok.end()) break;
			j = boost::lexical_cast<int>(*tokeniter);
		}
		else if (*tokeniter == "k")
		{
			++tokeniter;
      if ( tokeniter == tok.end()) break;
			k = boost::lexical_cast<int>(*tokeniter);
		}
	}

	return 0;
}
*/
