#include "PHG4SiliconTrackerSteppingAction.h"
#include "PHG4SiliconTrackerDetector.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Shower.h>

#include <g4main/PHG4TrackUserInfoV1.h>

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4TouchableHandle.hh>
#include <Geant4/G4TouchableHistory.hh>
#include <Geant4/G4ThreeVector.hh>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
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

#include <gsl/gsl_math.h>

#include <iostream>

//____________________________________________________________________________..
PHG4SiliconTrackerSteppingAction::PHG4SiliconTrackerSteppingAction(PHG4SiliconTrackerDetector* detector/*, const PHG4Parameters *parameters*/):
    detector_( detector ),
    hits_(NULL),
    absorberhits_(NULL),
    hit(NULL),
    // params(parameters),
    savehitcontainer(NULL),
    saveshower(NULL)
    // IsActive(params->get_int_param("active")),
    // IsBlackHole(params->get_int_param("blackhole"))
{
  std::cout << "PHG4SiliconTrackerSteppingAction created" << std::endl;
  verbosity = 0;
}

PHG4SiliconTrackerSteppingAction::~PHG4SiliconTrackerSteppingAction()
{
  // if the last hit was a zero energie deposit hit, it is just reset
  // and the memory is still allocated, so we need to delete it here
  // if the last hit was saved, hit is a NULL pointer which are
  // legal to delete (it results in a no operation)
  delete hit;
}

//____________________________________________________________________________..
bool PHG4SiliconTrackerSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{

  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  // get volume of the current step
  G4VPhysicalVolume* volume = touch->GetVolume();

  const int whichactive = detector_->IsInSiliconTracker(volume);

  if (!whichactive)
    return false;

  // set ladder index
  int sphxlayer     = 0;
  int inttlayer     = 0;
  int ladderz       = 0;
  int ladderphi     = 0;
  int strip_z_index = 0;
  int strip_y_index = 0;

  if (whichactive > 0) // silicon acrive sensor
    {
      if(verbosity > 1)
	{
	  std::cout << std::endl << "PHG4SilicoTrackerSteppingAction::UserSteppingAction for volume name (pre) " <<  touch->GetVolume()->GetName() 
		    << " volume->GetTranslation " << touch->GetVolume()->GetTranslation()
		    << " volume->GetCopyNo() " << volume->GetCopyNo() 
		    << std::endl;
	  G4TouchableHandle touch_post = aStep->GetPostStepPoint()->GetTouchableHandle();
	  G4VPhysicalVolume* volume_post = touch_post->GetVolume();
	  std::cout  << "PHG4SilicoTrackerSteppingAction::UserSteppingAction for volume name (post) " <<  touch_post->GetVolume()->GetName() 
		     << " volume->GetTranslation " << touch_post->GetVolume()->GetTranslation()
		     << " volume->GetCopyNo() " << volume_post->GetCopyNo() 
		     << std::endl;
	  std::cout << " IsFirstStepinVolume = " << aStep->IsFirstStepInVolume() << " IsLastStepInVolume = " << aStep->IsLastStepInVolume() << std::endl;
	}
      
      boost::char_separator<char> sep("_");
      boost::tokenizer<boost::char_separator<char> > tok(touch->GetVolume(2)->GetName(), sep);
      boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
      for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
	{
	  const int index = std::distance(tok.begin(), tokeniter);
	  // skip index=0 for preamble

	  if (index==1)
	    sphxlayer = boost::lexical_cast<int>(*tokeniter);

	  if (index==2)
	    inttlayer = boost::lexical_cast<int>(*tokeniter);

	  if (index==3)
	    ladderz   = boost::lexical_cast<int>(*tokeniter);

	  if (index==4)
	    ladderphi = boost::lexical_cast<int>(*tokeniter);
	}
      if (inttlayer < 0 || inttlayer >= 4)
	assert(!"PHG4SiliconTrackerSteppingAction: check INTT ladder layer.");

      // convert ladder type [0-3] to silicon sensor type [0-1]
      const int laddertype       = (ladderz==1 || ladderz==2) ? 0 : 1;
      const double strip_z       = (inttlayer==0) ?          detector_->arr_strip_z[0][laddertype] :          detector_->arr_strip_z[1][laddertype];
      const int nstrips_z_sensor = (inttlayer==0) ? detector_->arr_nstrips_z_sensor[0][laddertype] : detector_->arr_nstrips_z_sensor[1][laddertype];
      const int nstrips_phi_cell = detector_->arr_nstrips_phi_cell[inttlayer];
      const double strip_y       = detector_->arr_strip_y[inttlayer];

      // Find the strip y and z index values
      // This just regurgitates the values set in PHG4SiliconTrackerParameterization
      // when the G4PVParameterized was defined 
      G4ThreeVector strip_pos = volume->GetTranslation();

      strip_z_index = 0;
      for (int i=0; i<nstrips_z_sensor; ++i)
        {
          const double zmin = 2.*strip_z*(double)(i)   - strip_z*(double)nstrips_z_sensor;
          const double zmax = 2.*strip_z*(double)(i+1) - strip_z*(double)nstrips_z_sensor;
          if (strip_pos.z()/CLHEP::mm>zmin && strip_pos.z()/CLHEP::mm<=zmax)
            strip_z_index = i;
        }

      strip_y_index = 0;
      for (int i=0; i<2*nstrips_phi_cell; ++i)
        {
          const double ymin = 2.*strip_y*(double)(i)   - 2.*strip_y*(double)nstrips_phi_cell;
          const double ymax = 2.*strip_y*(double)(i+1) - 2.*strip_y*(double)nstrips_phi_cell;
          if (strip_pos.y()/CLHEP::mm>ymin && strip_pos.y()/CLHEP::mm<=ymax)
            {
	      strip_y_index = i;
	      //std::cout << " found strip y index = " << i << std::endl;
	      //std::cout << " i " << i << " ymin " << ymin << " ymax " << ymax << std::endl;
	    }
        }

      // The following is a hack to get around what seems to be a bug that causes the prestep to point to the wrong strip
      // The symptom of this is that this is the first and last step in the volume, AND the prestep and posttsep both have the same strip copy no.
      // There is no legitimate way for this to happen

      if( aStep->IsFirstStepInVolume() == 1 &&  aStep->IsLastStepInVolume() == 1)
	{
	  G4TouchableHandle touch_post = aStep->GetPostStepPoint()->GetTouchableHandle();
	  G4VPhysicalVolume* volume_post = touch_post->GetVolume();
	  
	  if(verbosity > 1)
	    {
	      std::cout << " ******** First and last step in this volume " << std::endl;
	      std::cout << "    volume_pre " <<  volume->GetName() << " volume_post " << volume_post->GetName() << std::endl;
	    }
	  if(volume->GetCopyNo() == volume_post->GetCopyNo())
	    {
	      if(verbosity > 1)
		std::cout <<        "     ************** Bogus! SAME COPY NUMBER FOR PRE AND POST VOLUMES MUST BE WRONG"  << std::endl;
	      
	      // we need a hack to replace the values above with the correct strip index values 
	      // the transform of the world coordinates into the sensor frame will work correctly, so we determine the strip indices from the hit position

	      G4StepPoint * prePoint = aStep->GetPreStepPoint();
	      G4ThreeVector preworldPos = prePoint->GetPosition();
	      G4ThreeVector strip_pos =  touch->GetHistory()->GetTransform(touch->GetHistory()->GetDepth() - 1).TransformPoint(preworldPos);
	  
	      strip_z_index = 0;
	      for (int i=0; i<nstrips_z_sensor; ++i)
		{
		  const double zmin = 2.*strip_z*(double)(i)   - strip_z*(double)nstrips_z_sensor;
		  const double zmax = 2.*strip_z*(double)(i+1) - strip_z*(double)nstrips_z_sensor;
		  if (strip_pos.z()/CLHEP::mm>zmin && strip_pos.z()/CLHEP::mm<=zmax)
		    strip_z_index = i;
		}
	      
	      strip_y_index = 0;
	      for (int i=0; i<2*nstrips_phi_cell; ++i)
		{
		  const double ymin = 2.*strip_y*(double)(i)   - 2.*strip_y*(double)nstrips_phi_cell;
		  const double ymax = 2.*strip_y*(double)(i+1) - 2.*strip_y*(double)nstrips_phi_cell;
		  if (strip_pos.y()/CLHEP::mm>ymin && strip_pos.y()/CLHEP::mm<=ymax)
		    {
		      strip_y_index = i;
		      //std::cout << " found strip y index = " << i << std::endl;
		      //std::cout << " i " << i << " ymin " << ymin << " ymax " << ymax << std::endl;
		    }
		}    
	      
	    }	  
	}
      /*
      std::cout << " sphxlayer " << sphxlayer
		<< " inttlayer " << inttlayer
		<< " ladderz " << ladderz
		<< " ladderphi " << ladderphi
		<< " strip_z_index " << strip_z_index 
		<< " strip_y_index " << strip_y_index
		<< " strip_pos.x " << strip_pos.x()
		<< " strip_pos.y " << strip_pos.y()
		<< " strip_pos.z " << strip_pos.z()
		<< std::endl;
      */
    }
  else // silicon inactive area, FPHX, stabe etc. as absorbers
    {
      sphxlayer     = -1;
      inttlayer     = -1;
      ladderz       = -1;
      ladderphi     = -1;
      strip_z_index = -1;
      strip_y_index = -1;
    }
  
  // collect energy and track length step by step
  G4double edep =  aStep->GetTotalEnergyDeposit()/CLHEP::GeV;
  G4double eion = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit())/CLHEP::GeV;
  const G4Track* aTrack = aStep->GetTrack();

  // if this block stops everything, just put all kinetic energy into edep
  if (/*IsBlackHole*/detector_->IsBlackHole())
    {
      edep = aTrack->GetKineticEnergy()/CLHEP::GeV;
      G4Track* killtrack = const_cast<G4Track *> (aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }

  // make sure we are in a volume
  if (/*IsActive*/detector_->IsActive())
    {
      bool geantino = false;

      // the check for the pdg code speeds things up, I do not want to make
      // an expensive string compare for every track when we know
      // geantino or chargedgeantino has pid=0
      if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 && aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos)
	geantino = true;

      G4StepPoint * prePoint = aStep->GetPreStepPoint();
      G4StepPoint * postPoint = aStep->GetPostStepPoint();
      if(verbosity > 1)
	std::cout << "prePoint step status = " << prePoint->GetStepStatus() << " postPoint step status = " << postPoint->GetStepStatus() << std::endl;
      switch (prePoint->GetStepStatus())
        {
        case fGeomBoundary:
        case fUndefined:

          // if previous hit was saved, hit pointer was set to NULL
          // and we have to make a new one
          if (! hit)
            {
              hit = new PHG4Hitv1();
            }

          hit->set_layer((unsigned int)sphxlayer);

          // set the index values needed to locate the sensor strip
          hit->set_strip_z_index(strip_z_index);
          hit->set_strip_y_index(strip_y_index);
          hit->set_ladder_z_index(ladderz);
          hit->set_ladder_phi_index(ladderphi);

          //here we set the entrance values in cm
          hit->set_x(0, prePoint->GetPosition().x()/CLHEP::cm);
          hit->set_y(0, prePoint->GetPosition().y()/CLHEP::cm);
          hit->set_z(0, prePoint->GetPosition().z()/CLHEP::cm);

          hit->set_px(0, prePoint->GetMomentum().x()/CLHEP::GeV);
          hit->set_py(0, prePoint->GetMomentum().y()/CLHEP::GeV);
          hit->set_pz(0, prePoint->GetMomentum().z()/CLHEP::GeV);

          // time in ns
          hit->set_t(0, prePoint->GetGlobalTime()/CLHEP::nanosecond);

          //set the track ID
          hit->set_trkid(aTrack->GetTrackID());

          //set the initial energy deposit
          hit->set_edep(0);
          hit->set_eion(0); // only implemented for v5 otherwise empty

          if (whichactive > 0) // return of IsInSiliconTracker, > 0 hit in si-strip, < 0 hit in absorber
            {
              // Now save the container we want to add this hit to
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

          break;

        default:
          break;
        }

      // here we just update the exit values, it will be overwritten
      // for every step until we leave the volume or the particle
      // ceases to exist
      hit->set_x( 1, postPoint->GetPosition().x()/CLHEP::cm);
      hit->set_y( 1, postPoint->GetPosition().y()/CLHEP::cm);
      hit->set_z( 1, postPoint->GetPosition().z()/CLHEP::cm);

      hit->set_px(1, postPoint->GetMomentum().x()/CLHEP::GeV);
      hit->set_py(1, postPoint->GetMomentum().y()/CLHEP::GeV);
      hit->set_pz(1, postPoint->GetMomentum().z()/CLHEP::GeV);

      hit->set_t( 1, postPoint->GetGlobalTime()/CLHEP::nanosecond);

      //sum up the energy to get total deposited
      hit->set_edep(hit->get_edep() + edep);
      hit->set_eion(hit->get_eion() + eion);

      if (geantino)
	{
	  hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression
	  hit->set_eion(-1);
	}

      if (edep > 0)
        if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
          if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
            pp->SetKeep(1); // we want to keep the track

      // if any of these conditions is true this is the last step in
      // this volume and we need to save the hit
      // postPoint->GetStepStatus() == fGeomBoundary: track leaves this volume
      // postPoint->GetStepStatus() == fWorldBoundary: track leaves this world
      // (not sure if this will ever be the case)
      // aTrack->GetTrackStatus() == fStopAndKill: track ends
      if (postPoint->GetStepStatus() == fGeomBoundary || postPoint->GetStepStatus() == fWorldBoundary|| aTrack->GetTrackStatus() == fStopAndKill)
        {
	  if(verbosity > 1)
	    std::cout << " postPoint step status changed, save hit and delete it" << std::endl;

          // save only hits with energy deposit (or -1 for geantino)
          if (hit->get_edep())
            {
              savehitcontainer->AddHit(sphxlayer, hit);
              if (saveshower)
                {
                  saveshower->add_g4hit_id(hits_->GetID(),hit->get_hit_id());
                }
	      if(verbosity > 1)
		hit->print();
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

      if(verbosity > 1)
	{
          G4StepPoint * prePoint = aStep->GetPreStepPoint();
          G4StepPoint * postPoint = aStep->GetPostStepPoint();
	  G4ThreeVector preworldPos = prePoint->GetPosition();
	  G4ThreeVector postworldPos = postPoint->GetPosition();

	  std::cout <<  " entry point world pos " <<  prePoint->GetPosition().x() << "  " <<  prePoint->GetPosition().y()  << "  " <<  prePoint->GetPosition().z() << std::endl; 
	  std::cout << " exit point world pos " << postPoint->GetPosition().x() << "  " <<  postPoint->GetPosition().y()  << "  " <<  postPoint->GetPosition().z() << std::endl; 

	  // The exit point transforms do not work here because the particle has already entered the next volume 
	  // - go back and find Jin's fix for this

	  // strip local pos
	  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(prePoint->GetTouchable());
	  G4ThreeVector prelocalPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(preworldPos);
	  std::cout << " entry point strip local pos: " << " is " << prelocalPos.x() << " " << prelocalPos.y() << " " << prelocalPos.z() << std::endl;
	  G4TouchableHistory* postTouchable = (G4TouchableHistory*)(postPoint->GetTouchable());
	  G4ThreeVector postlocalPos = postTouchable->GetHistory()->GetTopTransform().TransformPoint(postworldPos);
	  std::cout << " exit point strip local pos: " <<postlocalPos.x() << " " << postlocalPos.y() << " " << postlocalPos.z() << std::endl;

	  // sensor local pos
	  G4ThreeVector presensorLocalPos = theTouchable->GetHistory()->GetTransform(theTouchable->GetHistory()->GetDepth() - 1).TransformPoint(preworldPos);
	  std::cout << " entry point sensor local pos: " << presensorLocalPos.x() << " " << presensorLocalPos.y() << " " << presensorLocalPos.z() << std::endl;
	  G4ThreeVector postsensorLocalPos = postTouchable->GetHistory()->GetTransform(postTouchable->GetHistory()->GetDepth() - 1).TransformPoint(postworldPos);
	  std::cout << " exit point sensor local pos: " << postsensorLocalPos.x() << " " << postsensorLocalPos.y() << " " << postsensorLocalPos.z() << std::endl;
	}

      if (whichactive>0 && verbosity>0) // return of IsInSiliconTracker, > 0 hit in si-strip, < 0 hit in absorber
        {
          G4StepPoint * prePoint = aStep->GetPreStepPoint();
          G4StepPoint * postPoint = aStep->GetPostStepPoint();
          std::cout << "----- PHg4SiliconTrackerSteppingAction::UserSteppingAction - active volume = " << volume->GetName() << std::endl;
          std::cout << "       strip_z_index = " << strip_z_index << " strip_y_index = " << strip_y_index << std::endl;
          std::cout << "       prepoint x position "  << prePoint->GetPosition().x() /CLHEP::cm << std::endl;
          std::cout << "       prepoint y position "  << prePoint->GetPosition().y() /CLHEP::cm << std::endl;
          std::cout << "       prepoint z position "  << prePoint->GetPosition().z() /CLHEP::cm << std::endl;
          std::cout << "       postpoint x position " << postPoint->GetPosition().x()/CLHEP::cm << std::endl;
          std::cout << "       postpoint y position " << postPoint->GetPosition().y()/CLHEP::cm << std::endl;
          std::cout << "       postpoint z position " << postPoint->GetPosition().z()/CLHEP::cm << std::endl;
          std::cout << "       edep " << edep << std::endl;
          std::cout << "       eion " << eion << std::endl;
        }

      return true;
    }
  else
    {
      return false;
    }
}

//____________________________________________________________________________..
void PHG4SiliconTrackerSteppingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  const std::string detectorname = (detector_->SuperDetector() != "NONE") ? detector_->SuperDetector() : detector_->GetName();
  const std::string hitnodename      = "G4HIT_"          + detectorname;
  const std::string absorbernodename = "G4HIT_ABSORBER_" + detectorname;

  //now look for the map and grab a pointer to it.
  hits_         =  findNode::getClass<PHG4HitContainer>(topNode , hitnodename.c_str());
  absorberhits_ =  findNode::getClass<PHG4HitContainer>(topNode , absorbernodename.c_str());

  // if we do not find the node it's messed up.
  if (!hits_)
    std::cout << "PHG4SiliconTrackerSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;

  if (!absorberhits_ && verbosity>1)
    std::cout << "PHG4SiliconTrackerSteppingAction::SetTopNode - unable to find " << absorbernodename << std::endl;
}
