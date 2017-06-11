#include <g4main/PHG4TrackUserInfoV1.h>
#include <phool/getClass.h>
#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <TMath.h>

#include <iomanip>
#include <iostream>

#include "TPCbase/TPCHit.h"
#include "TPCbase/TPCHitsContainer.h"
#include "TPCbase/TPCConstants.h"
#include "TPCSteppingAction.h"
#include "TPCDetector.h"

#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>

using namespace std;
using namespace TPCConstants;
//____________________________________________________________________________..
TPCSteppingAction::TPCSteppingAction(TPCDetector* detector):
  fVerbosity(0),
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
  if(fPHHits) SteppingActionPHG4(aStep);
  return SteppingActionTPC(aStep);
}
//____________________________________________________________________________..
bool TPCSteppingAction::SteppingActionTPC( const G4Step *aStep) {
  const G4Track* aTrack = aStep->GetTrack();
  G4ParticleDefinition *partDef = aTrack->GetDefinition();
  if(partDef->GetPDGCharge()==0) return false; // skip neutral particles
  //if(partDef->GetPDGEncoding()==11) return false; // skip electrons (wanted to skip only deltas, but ...)

  if(fVerbosity>10) {
    std::cout << " Track from: " << partDef->GetParticleName() << " " << partDef->GetPDGEncoding() << std::endl;
    std::cout << " IsGoodForTracking " << (aTrack->IsGoodForTracking()?"yes":"no ") << std::endl;
    std::cout << " length " << aStep->GetStepLength()/cm << std::endl;
  }
  fHit->SetTrack( aTrack->GetTrackID() );
  //--- G4Step ==> TPCHit
  G4double edep = aStep->GetTotalEnergyDeposit() - aStep->GetNonIonizingEnergyDeposit();
  fHit->SetDEnergy(edep/keV);
  G4StepPoint *prePoint = aStep->GetPreStepPoint();
  Float_t x0 = prePoint->GetPosition().x()/cm;
  Float_t y0 = prePoint->GetPosition().y()/cm;
  Float_t z0 = prePoint->GetPosition().z()/cm;
  //fHit->set_px( 0, prePoint->GetMomentum().x() / GeV );
  //fHit->set_py( 0, prePoint->GetMomentum().y() / GeV );
  //fHit->set_pz( 0, prePoint->GetMomentum().z() / GeV );
  //fHit->set_t( 0, prePoint->GetGlobalTime() / nanosecond );
  G4StepPoint *postPoint = aStep->GetPostStepPoint();
  Float_t x1 = postPoint->GetPosition().x()/cm;
  Float_t y1 = postPoint->GetPosition().y()/cm;
  Float_t z1 = postPoint->GetPosition().z()/cm;
  Float_t dx = x1-x0;
  Float_t dy = y1-y0;
  Float_t dz = z1-z0;
  fHit->SetL( TMath::Sqrt( dx*dx + dy*dy + dz*dz ) );
  Float_t x = 0.5*(x0+x1);
  Float_t y = 0.5*(y0+y1);
  Float_t z = 0.5*(z0+z1);
  fHit->SetZ( z );
  fHit->SetR( TMath::Sqrt( x*x + y*y ) );
  fHit->SetPhi( TMath::ATan2(-y,-x)+TMath::Pi() );
  if(fVerbosity>10) {
    std::cout << " R" << fHit->GetR() <<
      " P" << fHit->GetPhi() <<  " Z" << fHit->GetZ() << " L" << fHit->GetL() <<
      " E" << fHit->GetDEnergy() << " TID" << fHit->GetTrack() << std::endl;
    std::cout << " Rextend[" << TMath::Sqrt(x0*x0+y0*y0) << "-" << TMath::Sqrt(x1*x1+y1*y1) << "]" << std::endl;
  }
  //---
  if(fHit->GetDEnergy()>0) fHits->Add(fHit); // container does not own fHit!
  return true;
}
//____________________________________________________________________________..
void TPCSteppingAction::SteppingActionPHG4( const G4Step *aStep) {
  //saving the whole step for compatibility to existing external tools
  const G4Track* aTrack = aStep->GetTrack();
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
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
  float rads = hit->get_avg_x()*hit->get_avg_x();
  rads += hit->get_avg_y()*hit->get_avg_y();
  rads = TMath::Sqrt( rads );
  unsigned int layers = kNPadRowsPerModule * kNSections;
  float steplyr = (kGasOuterRadius-kGasOuterRadius)/int(layers);
  unsigned int lyr = (rads-kGasInnerRadius)/steplyr;
  hit->set_layer( lyr );
  if(lyr>layers) lyr = layers;
  fPHHits->AddHit(lyr, hit); // container own hit from now on.
}
//=====
void TPCSteppingAction::SetInterfacePointers( PHCompositeNode* node ) {
  fHits =  findNode::getClass<TPCHitsContainer>( node , "TPCHits" );
  fPHHits =  findNode::getClass<PHG4HitContainer>( node , "G4HIT_TPC" );
}
