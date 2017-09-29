#include <g4main/PHG4TrackUserInfoV1.h>
#include <phool/getClass.h>
#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <TMath.h>

#include <iomanip>
#include <iostream>

#include <TPCHit.h>
#include <TPCHitsContainer.h>
#include <TPCConstants.h>
#include "TPCSteppingAction.h"
#include "TPCDetector.h"

#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>

using namespace std;
using namespace TPCConstants;
//____________________________________________________________________________..
TPCSteppingAction::TPCSteppingAction(TPCDetector* detector):
  fVerbosity(0),
  fSkipNeutral(false),
  fDetector( detector ),
  fHits(NULL),
  fHit(new TPCHit()),
  fPHHits(NULL)
{
}

TPCSteppingAction::~TPCSteppingAction() {
  delete fHit;
}

//____________________________________________________________________________..
bool TPCSteppingAction::UserSteppingAction( const G4Step* aStep, bool )
{
  G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if(!fDetector->IsThisActive(volume)) return false;
  if(fVerbosity>10) {
    std::cout << "TPCSteppingAction::UserSteppingAction" << std::endl;
    std::cout << volume->GetName() << std::endl;
  }
  return SteppingActionTPC(aStep);
}
//____________________________________________________________________________..
bool TPCSteppingAction::SteppingActionTPC( const G4Step *aStep) {
  const G4Track* aTrack = aStep->GetTrack();
  G4ParticleDefinition *partDef = aTrack->GetDefinition();
  if(fSkipNeutral) if(partDef->GetPDGCharge()==0) return false; // skip neutral particles
  if(fVerbosity>10) {
    std::cout << " Track from: " << partDef->GetParticleName() << " " << partDef->GetPDGEncoding() << std::endl;
    std::cout << " IsGoodForTracking " << (aTrack->IsGoodForTracking()?"yes":"no ") << std::endl;
    std::cout << " length " << aStep->GetStepLength()/cm << std::endl;
  }

  // reads from G4step
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();

  // fills data container
  PHG4Hitv1 *hit = new PHG4Hitv1();
  hit->set_x(0, prePoint->GetPosition().x() / cm);
  hit->set_y(0, prePoint->GetPosition().y() / cm);
  hit->set_z(0, prePoint->GetPosition().z() / cm);
  hit->set_px(0, prePoint->GetMomentum().x() / GeV);
  hit->set_py(0, prePoint->GetMomentum().y() / GeV);
  hit->set_pz(0, prePoint->GetMomentum().z() / GeV);
  hit->set_t(0, prePoint->GetGlobalTime() / nanosecond);
  hit->set_trkid(aTrack->GetTrackID());
  hit->set_x(1, postPoint->GetPosition().x() / cm);
  hit->set_y(1, postPoint->GetPosition().y() / cm);
  hit->set_z(1, postPoint->GetPosition().z() / cm);
  hit->set_px(1, postPoint->GetMomentum().x() / GeV);
  hit->set_py(1, postPoint->GetMomentum().y() / GeV);
  hit->set_pz(1, postPoint->GetMomentum().z() / GeV);
  hit->set_t(1, postPoint->GetGlobalTime() / nanosecond);
  hit->set_edep( aStep->GetTotalEnergyDeposit()/GeV );
  G4double edep = aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit();
  hit->set_eion(edep/GeV);
  float rads = hit->get_avg_x()*hit->get_avg_x();
  rads += hit->get_avg_y()*hit->get_avg_y();
  rads = TMath::Sqrt( rads );
  unsigned int layers = kNPadRowsPerModule * kNSections;
  float steplyr = (kGasOuterRadius-kGasOuterRadius)/int(layers);
  unsigned int lyr = (rads-kGasInnerRadius)/steplyr;
  hit->set_layer( lyr );
  if(lyr>layers) lyr = layers;
  fPHHits->AddHit(lyr, hit); // container own hit from now on.
  // add hit to sorter
  fHit->Assign( hit );
  if(fVerbosity>10) {
    std::cout << " X" << fHit->GetX() <<
      " Y" << fHit->GetY() <<  " Z" << fHit->GetZ() << " L" << fHit->GetL() <<
      " E" << fHit->GetDEnergy() << " TID" << fHit->GetTrack() << std::endl;
  }
  if(fHit->GetDEnergy()>0) fHits->Add(fHit); // container does not own fHit!
  return true;
}
//=====
void TPCSteppingAction::SetInterfacePointers( PHCompositeNode* node ) {
  fPHHits =  findNode::getClass<PHG4HitContainer>( node , "G4HIT_SVTX" );
  fHits =  findNode::getClass<TPCHitsContainer>( node , "TPCHits" );
}
