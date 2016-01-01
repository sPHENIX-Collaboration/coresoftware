#include "PHG4SiliconTrackerSteppingAction.h"
#include "PHG4SiliconTrackerDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

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
PHG4SiliconTrackerSteppingAction::PHG4SiliconTrackerSteppingAction( PHG4SiliconTrackerDetector* detector ):
  PHG4SteppingAction(NULL),
  detector_( detector ),
  hits_(NULL),
  absorberhits_(NULL),
  hit(NULL)
{
  cout << "PHG4SiliconTrackerSteppingAction created" << endl;
}

//____________________________________________________________________________..
bool PHG4SiliconTrackerSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* strip_volume = touch->GetVolume();

  // PHG4SiliconTrackerDetector_->IsInSiliconTracker(volume)
  // returns 
  //  0 if outside of SiliconTracker
  //  1 if inside sensor

  int whichactive = detector_->IsInSiliconTracker(strip_volume);

  if (!whichactive)
    {
      return false;
    }

  /*
  // Look at the ancestry of this volume
  G4VPhysicalVolume* v = touch->GetVolume();
  cout << "volume  is " << v->GetName() << endl;
  G4VPhysicalVolume* v1 = touch->GetVolume(1);
  cout << "volume 1 up is " << v1->GetName() << endl;
  G4VPhysicalVolume* v2 = touch->GetVolume(2);
  cout << "volume 2 up is " << v2->GetName() << endl;
  G4VPhysicalVolume* v3 = touch->GetVolume(3);
  cout << "volume 3 up is " << v3->GetName() << endl;
  */

  //===========================
  // we want the location of the hit strip
  //===========================

  // The location of the strip in its sensor (one sensor per ladder segment) can be obstained from the strip name
  // eg.  "av_2_impr_1_layer_0_sensor_strip_9_170_pv_15074" is the strip at icolumn = 9 and istrip = 170 

  int strip_column = -1;
  int strip_index = -1;

  // FYI: doing string compares inside a stepping action sounds like a recipe
  // for failure inside a heavy ion event... we'll wait and see how badly
  // this profiles. -MPM
  
  boost::char_separator<char> sep("_");
  boost::tokenizer<boost::char_separator<char> > tok(strip_volume->GetName(), sep);
  boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
  for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
    {
      if (*tokeniter == "address")
	{
	  // the next field is the strip column
	  ++tokeniter;
	  strip_column = boost::lexical_cast<int>(*tokeniter);
	  
	  // now the strip index in the column
	  ++tokeniter;
	  strip_index = boost::lexical_cast<int>(*tokeniter);
	}
    }
  
  // The ladder segment is the grandmother of the strip
  // ladder segment name contains the z and phi index of the segment
  // eg. "ladder_segment_6_17" is the ladder at iz = 6 and iphi = 17

  G4VPhysicalVolume *ladder_segment = touch->GetVolume(2);
  
  int segment_z_index = -1;
  int segment_phi_index = -1;
  boost::tokenizer<boost::char_separator<char> > tok1(ladder_segment->GetName(), sep);
  for (tokeniter = tok1.begin(); tokeniter != tok1.end(); ++tokeniter)
    {
      if(*tokeniter == "segment")
	{	  
	  // the next field is the z index
	  ++tokeniter;
	  segment_z_index = boost::lexical_cast<int>(*tokeniter);
	  
	  // now the phi index
	  ++tokeniter;
	  segment_phi_index = boost::lexical_cast<int>(*tokeniter);
	}
    }      
  
  // Now we want to collect information about the hit
  
  G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
  const G4Track* aTrack = aStep->GetTrack();
  
  // if this cylinder stops everything, just put all kinetic energy into edep
  if (detector_->IsBlackHole())
    {
      edep = aTrack->GetKineticEnergy()/GeV;
      G4Track* killtrack = const_cast<G4Track *> (aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }
  
  int layer_id = detector_->get_Layer();
  //cout << endl << "  In UserSteppingAction for layer " << layer_id << endl;
  
  // test if we are active
  if ( detector_->IsActive() )
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
      //        cout << "time prepoint: " << prePoint->GetGlobalTime()/ns << endl;
      //        cout << "time postpoint: " << postPoint->GetGlobalTime()/ns << endl;
      //        cout << "kinetic energy: " <<  aTrack->GetKineticEnergy()/GeV << endl;
      //       G4ParticleDefinition* def = aTrack->GetDefinition();
      //       cout << "Particle: " << def->GetParticleName() << endl;
      switch (prePoint->GetStepStatus())
	{
	case fGeomBoundary:
	case fUndefined:
	  
	  hit = new PHG4Hitv1();
	  
	  hit->set_layer((unsigned int)layer_id);
	  
	  // set the index values needed to locate the sensor strip
	  hit->set_strip_z_index(strip_column);
	  hit->set_strip_y_index(strip_index);
	  hit->set_ladder_z_index(segment_z_index);
	  hit->set_ladder_phi_index(segment_phi_index);
	  
	  
	  //here we set the entrance values in cm
	  hit->set_x( 0, prePoint->GetPosition().x() / cm );
	  hit->set_y( 0, prePoint->GetPosition().y() / cm );
	  hit->set_z( 0, prePoint->GetPosition().z() / cm );
	  
	  hit->set_px( 0, prePoint->GetMomentum().x() / GeV );
	  hit->set_py( 0, prePoint->GetMomentum().y() / GeV );
	  hit->set_pz( 0, prePoint->GetMomentum().z() / GeV );
	  
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
	  
	  // Now add the hit
	  //	  hit->print();
	  hits_->AddHit(layer_id, hit);
	  
	  {
	    if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
	      {
		if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
		  {
		    pp->GetShower()->add_g4hit_id(hits_->GetID(),hit->get_hit_id());
		  }
	      }
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
      
      hit->set_px(1, postPoint->GetMomentum().x() / GeV );
      hit->set_py(1, postPoint->GetMomentum().y() / GeV );
      hit->set_pz(1, postPoint->GetMomentum().z() / GeV );
      
      hit->set_t( 1, postPoint->GetGlobalTime() / nanosecond );
      //sum up the energy to get total deposited
      hit->set_edep(hit->get_edep() + edep);
      if (geantino)
	{
	  hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
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
      
      if(verbosity>0)
	{  
	  G4StepPoint * prePoint = aStep->GetPreStepPoint();
	  G4StepPoint * postPoint = aStep->GetPostStepPoint();
	  cout << "----- PHg4SiliconTrackerSteppingAction::UserSteppingAction - active volume = " << strip_volume->GetName() << endl;
	  cout << "       strip_column = " << strip_column << " strip_index = " << strip_index << endl;
	  cout << "       ladder segment name = " << ladder_segment->GetName() << endl;
	  cout << "       segment_z_index = " << segment_z_index << " segment_phi_index = " << segment_phi_index << endl;
	  cout << "       prepoint x position " << prePoint->GetPosition().x() / cm << endl;
	  cout << "       prepoint y position " << prePoint->GetPosition().y() / cm << endl;
	  cout << "       prepoint z position " << prePoint->GetPosition().z() / cm << endl;
	  cout << "       postpoint x position " << postPoint->GetPosition().x() / cm << endl;
	  cout << "       postpoint y position " << postPoint->GetPosition().y() / cm << endl;
	  cout << "       postpoint z position " << postPoint->GetPosition().z() / cm << endl;
	  cout << "       edep " << edep  << endl;
	}

      // way to chatty
      //cout << "  stepping action found hit:" << endl;
      //hit->identify(); 
      //cout << endl;
      
      // return true to indicate the hit was used
      return true;
    }
  else
    {
      return false;
    }
  
} 

//____________________________________________________________________________..
void PHG4SiliconTrackerSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
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
      std::cout << "PHG4SiliconTrackerSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
    }
  if ( ! absorberhits_)
    {
      if (verbosity > 0)
	{
	  cout << "PHG4SiliconTrackerSteppingAction::SetTopNode - unable to find " << absorbernodename << endl;
	}
    }
}
