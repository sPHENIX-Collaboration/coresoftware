#include "PHG4CrystalCalorimeterSteppingAction.h"
#include "PHG4CrystalCalorimeterDetector.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv8.h>

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
  hit(NULL)
{

}


//____________________________________________________________________________..
bool PHG4CrystalCalorimeterSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{
  G4TouchableHandle touch = aStep->GetPreStepPoint()->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch->GetVolume();

  if (! detector_->IsInCrystalCalorimeter(volume) )
    {
      return false;
    }
  else
    {
      /* Find indizes of crystal containing this step */
      int layer_id = detector_->get_Layer();
      int idx_j = -1;
      int idx_k = -1;
      int idx_l = -1;

      boost::char_separator<char> sep("_");
      boost::tokenizer<boost::char_separator<char> > tok(volume->GetName(), sep);
      boost::tokenizer<boost::char_separator<char> >::const_iterator tokeniter;
      for (tokeniter = tok.begin(); tokeniter != tok.end(); ++tokeniter)
        {
          if (*tokeniter == "j")
            {
              ++tokeniter;
              idx_j = boost::lexical_cast<int>(*tokeniter);
            }
          else if (*tokeniter == "k")
            {
              ++tokeniter;
              idx_k = boost::lexical_cast<int>(*tokeniter);
            }
        }

      /* Get pointer to associated Geant4 track */
      const G4Track* aTrack = aStep->GetTrack();

      /* Check if particle is 'geantino' */
      bool geantino = false;
      if (aTrack->GetParticleDefinition()->GetPDGEncoding() == 0 &&
          aTrack->GetParticleDefinition()->GetParticleName().find("geantino") != string::npos)
        {
          geantino = true;
        }

      /* If inside active crystal and particle is NOT a geantino, add its energy deposited to output collection */
      if ( !geantino && detector_->IsActive() )
        {
          /* Get energy deposited by this step */
          G4double edep = aStep->GetTotalEnergyDeposit() / GeV;
          //G4double edep = (aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit()) / GeV;

          /* Get Geant4 pre- and post-step points */
          G4StepPoint * prePoint = aStep->GetPreStepPoint();
          G4StepPoint * postPoint = aStep->GetPostStepPoint();
          switch (prePoint->GetStepStatus())
            {
            case fGeomBoundary:
            case fUndefined:
              hit = new PHG4Hitv8();

              /* Set hit location (tower index) */
              hit->set_index_j(idx_j);
              hit->set_index_k(idx_k);
              hit->set_index_l(idx_l);

              /* Set hit location (space point) */
              hit->set_x( 0, prePoint->GetPosition().x() / cm);
              hit->set_y( 0, prePoint->GetPosition().y() / cm );
              hit->set_z( 0, prePoint->GetPosition().z() / cm );

              hit->set_x( 1, postPoint->GetPosition().x() / cm );
              hit->set_y( 1, postPoint->GetPosition().y() / cm );
              hit->set_z( 1, postPoint->GetPosition().z() / cm );

              /* Set hit time */
              hit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );
              hit->set_t( 1, postPoint->GetGlobalTime() / nanosecond );

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

                /* set hit energy */
                hit->set_edep( edep );
              }

              /* Now add the hit to the hit collection */
              hits_->AddHit(layer_id, hit);
              break;
            default:
              break;
            }

          return true;

        }
      else
        {
          return false;
        }
    }

  return true;

}


//____________________________________________________________________________..
void PHG4CrystalCalorimeterSteppingAction::SetInterfacePointers( PHCompositeNode* topNode )
{
  string hitnodename;
  hitnodename = "G4HIT_" + detector_->GetName();

  hits_ =  findNode::getClass<PHG4HitContainer>( topNode , hitnodename.c_str() );

  /* Check if we found a node */
  if ( ! hits_ )
    {
      std::cout << "PHG4CrystalCalorimeterSteppingAction::SetTopNode - unable to find " << hitnodename << std::endl;
    }
}
