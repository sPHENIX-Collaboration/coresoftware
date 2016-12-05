#include "PHG4SiliconTrackerSteppingAction.h"
#include "PHG4SiliconTrackerDetector.h"

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

#include <TMath.h>
#include <TVector3.h>
#include <TRotation.h>

#include <iostream>

//____________________________________________________________________________..
PHG4SiliconTrackerSteppingAction::PHG4SiliconTrackerSteppingAction(PHG4SiliconTrackerDetector* detector):
    detector_( detector ),
    hits_(NULL),
    absorberhits_(NULL),
    hit(NULL)
{
  std::cout << "PHG4SiliconTrackerSteppingAction created" << std::endl;
}

//____________________________________________________________________________..
bool PHG4SiliconTrackerSteppingAction::UserSteppingAction(const G4Step* aStep, bool)
{
  G4StepPoint            *preStep = aStep->GetPreStepPoint();
  G4TouchableHistory   *touchable = (G4TouchableHistory*)(preStep->GetTouchable());
  G4VPhysicalVolume *strip_volume = touchable->GetVolume();

  if (!detector_->IsInSiliconTracker(strip_volume))
    return false;

  // Now we want to collect information about the hit
  const G4Track* aTrack = aStep->GetTrack();
  G4double edep = aStep->GetTotalEnergyDeposit()/CLHEP::GeV;

  // if this cylinder stops everything, just put all kinetic energy into edep
  if (detector_->IsBlackHole())
    {
      edep = aTrack->GetKineticEnergy()/CLHEP::GeV;
      G4Track* killtrack = const_cast<G4Track *>(aTrack);
      killtrack->SetTrackStatus(fStopAndKill);
    }
  if (edep==0.)
    return false;

  const G4double nonionedep = aStep->GetNonIonizingEnergyDeposit();

  // set ladder index
  int sphxlayer = 0;
  int inttlayer = 0;
  int ladderz   = 0;
  int ladderphi = 0;

  boost::char_separator<char> sep("_");
  boost::tokenizer<boost::char_separator<char> > tok(touchable->GetVolume(2)->GetName(), sep);
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
  if (inttlayer < 0 || inttlayer >= 5)
    assert(!"PHG4SiliconTrackerSteppingAction: check INTT ladder layer.");

  // convert ladder type [0-3] to silicon sensor type [0-1]
  const int laddertype       = (ladderz==1 || ladderz==2) ? 0 : 1;
  const double strip_z       = (inttlayer==0) ?          detector_->arr_strip_z[0][laddertype] :          detector_->arr_strip_z[1][laddertype];
  const int nstrips_z_sensor = (inttlayer==0) ? detector_->arr_nstrips_z_sensor[0][laddertype] : detector_->arr_nstrips_z_sensor[1][laddertype];
  const int nstrips_phi_cell = detector_->arr_nstrips_phi_cell[inttlayer];
  const double strip_y       = detector_->arr_strip_y[inttlayer];

  G4ThreeVector strip_pos = strip_volume->GetTranslation();
  const double epsz = TMath::Sign(TMath::Power(10., -10.), strip_pos.z());
  const double epsy = TMath::Sign(TMath::Power(10., -10.), strip_pos.y());
  const int strip_z_index = (nstrips_z_sensor%2==0) ? int((strip_pos.z()/CLHEP::mm-strip_z+epsz)/(2.*strip_z)) + int(nstrips_z_sensor/2) : int((strip_pos.z()/CLHEP::mm)/(2.*strip_z)) + int(nstrips_z_sensor/2);
  const int strip_y_index = (nstrips_phi_cell%2==0) ? int((strip_pos.y()/CLHEP::mm-strip_y+epsy)/(2.*strip_y)) +     nstrips_phi_cell    : int((strip_pos.y()/CLHEP::mm)/(2.*strip_y)) +     nstrips_phi_cell;

  // test if we are active
  if (detector_->IsActive())
    {
      bool geantino = false;
      // the check for the pdg code speeds things up, I do not want to make
      // an expensive string compare for every track when we know
      // geantino or chargedgeantino has pid=0
      if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 && aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != std::string::npos)
        geantino = true;

      //G4StepPoint * preStep = aStep->GetPreStepPoint();
      switch (preStep->GetStepStatus())
        {
        case fGeomBoundary:

        case fUndefined:
          // create PHG4Hitv1's object
          hit = new PHG4Hitv1();

          hit->set_layer((unsigned int)sphxlayer);

          // set the index values needed to locate the sensor strip
          hit->set_strip_z_index(strip_z_index);
          hit->set_strip_y_index(strip_y_index);
          hit->set_ladder_z_index(ladderz);
          hit->set_ladder_phi_index(ladderphi);

          //here we set the entrance values in cm
          hit->set_x(0, preStep->GetPosition().x()/CLHEP::cm);
          hit->set_y(0, preStep->GetPosition().y()/CLHEP::cm);
          hit->set_z(0, preStep->GetPosition().z()/CLHEP::cm);

          hit->set_px(0, preStep->GetMomentum().x()/CLHEP::GeV);
          hit->set_py(0, preStep->GetMomentum().y()/CLHEP::GeV);
          hit->set_pz(0, preStep->GetMomentum().z()/CLHEP::GeV);

          // time in ns
          hit->set_t(0, preStep->GetGlobalTime()/CLHEP::nanosecond);

          //set the track ID
          {
            hit->set_trkid(aTrack->GetTrackID());
            if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
              if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
                {
                  hit->set_trkid(pp->GetUserTrackId());
                  hit->set_shower_id(pp->GetShower()->get_id());
                }
          }

          //set the initial energy deposit
          hit->set_edep(0);

          // Now add the hit
          /*hit->print();*/
          hits_->AddHit(sphxlayer, hit);

          {
            if (G4VUserTrackInformation* p = aTrack->GetUserInformation())
              if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
                pp->GetShower()->add_g4hit_id(hits_->GetID(),hit->get_hit_id());
          }

          break;

        default:
          break;
        }

      // here we just update the exit values, it will be overwritten
      // for every step until we leave the volume or the particle
      // ceases to exist
      G4StepPoint * postStep = aStep->GetPostStepPoint();
      hit->set_x( 1, postStep->GetPosition().x()/CLHEP::cm);
      hit->set_y( 1, postStep->GetPosition().y()/CLHEP::cm);
      hit->set_z( 1, postStep->GetPosition().z()/CLHEP::cm);

      hit->set_px(1, postStep->GetMomentum().x()/CLHEP::GeV);
      hit->set_py(1, postStep->GetMomentum().y()/CLHEP::GeV);
      hit->set_pz(1, postStep->GetMomentum().z()/CLHEP::GeV);

      hit->set_t( 1, postStep->GetGlobalTime()/CLHEP::nanosecond);

      //sum up the energy to get total deposited
      hit->set_edep(hit->get_edep() + edep);
      if (geantino)
        hit->set_edep(-1); // only energy=0 g4hits get dropped, this way geantinos survive the g4hit compression

      if (edep > 0)
        if ( G4VUserTrackInformation* p = aTrack->GetUserInformation() )
          if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) )
            pp->SetKeep(1); // we want to keep the track

      //const int verbosity = 0;
      if (verbosity>0)
        {
          G4StepPoint * preStep = aStep->GetPreStepPoint();
          G4StepPoint * postStep = aStep->GetPostStepPoint();
          std::cout << "----- PHg4SiliconTrackerSteppingAction::UserSteppingAction - active volume = " << strip_volume->GetName() << std::endl;
          std::cout << "       strip_z_index = " << strip_z_index << " strip_y_index = " << strip_y_index << std::endl;
          std::cout << "       prepoint x position "  << preStep->GetPosition().x() /CLHEP::cm << std::endl;
          std::cout << "       prepoint y position "  << preStep->GetPosition().y() /CLHEP::cm << std::endl;
          std::cout << "       prepoint z position "  << preStep->GetPosition().z() /CLHEP::cm << std::endl;
          std::cout << "       postpoint x position " << postStep->GetPosition().x()/CLHEP::cm << std::endl;
          std::cout << "       postpoint y position " << postStep->GetPosition().y()/CLHEP::cm << std::endl;
          std::cout << "       postpoint z position " << postStep->GetPosition().z()/CLHEP::cm << std::endl;
          std::cout << "       edep " << edep << std::endl;
          std::cout << "       nonionedep " << nonionedep << std::endl;
        }

      return true;
    }

  return false;
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

  if (!absorberhits_ && verbosity>0)
    std::cout << "PHG4SiliconTrackerSteppingAction::SetTopNode - unable to find " << absorbernodename << std::endl;
}
