#include "PHG4TruthTrackingAction.h"

#include "PHG4TruthEventAction.h"
#include "PHG4TrackUserInfoV1.h"
#include "PHG4UserPrimaryParticleInformation.h"

#include "PHG4Particlev2.h"
#include "PHG4VtxPointv1.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4Showerv1.h"

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackingManager.hh>
#include <Geant4/G4TrackVector.hh>
#include <Geant4/G4PrimaryParticle.hh>
#include <Geant4/G4VUserPrimaryParticleInformation.hh>

using namespace std;

const int VERBOSE = 0;

//________________________________________________________
PHG4TruthTrackingAction::PHG4TruthTrackingAction( PHG4TruthEventAction* eventAction ) :
  eventAction_( eventAction ), 
  truthInfoList_( NULL )
{}

void PHG4TruthTrackingAction::PreUserTrackingAction( const G4Track* track) {
   
  int trackid = 0;
  if (track->GetParentID()) {
    // secondaries get negative user ids and increment downward between geant subevents
    trackid = truthInfoList_->mintrkindex() - 1;
  } else {
    // primaries get positive user ids and increment upward between geant subevents
    trackid = truthInfoList_->maxtrkindex() + 1;
  }

  // add the user id to the geant4 user info
  PHG4TrackUserInfo::SetUserTrackId(const_cast<G4Track *> (track), trackid);
  
  // determine the momentum vector
  G4ParticleDefinition* def = track->GetDefinition();
  int pdgid = def->GetPDGEncoding();
  double m = def->GetPDGMass();
  double ke = track->GetVertexKineticEnergy();
  double ptot = sqrt(ke * ke + 2.0 * m * ke);

  G4ThreeVector pdir = track->GetVertexMomentumDirection();
  pdir *= ptot;

  // create a new particle -----------------------------------------------------
  PHG4Particlev2* ti = new PHG4Particlev2;
  ti->set_px(pdir[0] / GeV);
  ti->set_py(pdir[1] / GeV);
  ti->set_pz(pdir[2] / GeV);
  ti->set_track_id( trackid );

  ti->set_parent_id(track->GetParentID());
  if ( PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()) ) {
    ti->set_parent_id( p->GetUserParentId() );
  }

  ti->set_primary_id(trackid);
  if ( PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()) ) {
    if (track->GetParentID()) {
      ti->set_primary_id( p->GetUserPrimaryId() );
    } else {
      PHG4TrackUserInfo::SetUserPrimaryId(const_cast<G4Track *> (track), trackid);
      ti->set_primary_id(trackid);
    }
  }
  
  ti->set_pid( pdgid );
  ti->set_name(def->GetParticleName());
  ti->set_e(track->GetTotalEnergy() / GeV);

  if (!track->GetParentID()) {
    // primary track - propagate the barcode information
    PHG4UserPrimaryParticleInformation* userdata = static_cast<PHG4UserPrimaryParticleInformation*> (track->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
    if(userdata) ti->set_barcode( userdata->get_user_barcode() );
  }


  // create a new vertex object ------------------------------------------------
  G4ThreeVector v = track->GetVertexPosition();
  map<G4ThreeVector, int>::const_iterator viter = VertexMap.find(v);
  int vtxindex = 0;
  if (viter != VertexMap.end()) {
    vtxindex = viter->second;
  } else {

    vtxindex = truthInfoList_->maxvtxindex() + 1;
    if (track->GetParentID()) {
      vtxindex = truthInfoList_->minvtxindex() - 1;
    }

    VertexMap[v] = vtxindex;
    PHG4VtxPointv1 *vtxpt = new PHG4VtxPointv1(v[0] / cm,
					       v[1] / cm,
					       v[2] / cm,
					       track->GetGlobalTime() / ns);
    // insert new vertex into the output
    truthInfoList_->AddVertex(vtxindex, vtxpt);
  }

  ti->set_vtx_id(vtxindex);

  // insert particle into the output
  truthInfoList_->AddParticle(trackid, ti);


  // create or add to a new shower object --------------------------------------
  if (!track->GetParentID()) {
    PHG4Showerv1* shower = new PHG4Showerv1();
    PHG4TrackUserInfo::SetShower(const_cast<G4Track *> (track), shower);
    truthInfoList_->AddShower(trackid, shower);
    shower->set_id(trackid); // fyi, secondary showers may not share these ids
    shower->set_parent_particle_id(trackid);
    shower->set_parent_shower_id(0);
  } else {
    // get shower
    if ( G4VUserTrackInformation* p = track->GetUserInformation() ) {
      if ( PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p) ) {
	if (pp->GetShower()) {
	  pp->GetShower()->add_g4particle_id(trackid);
	  pp->GetShower()->add_g4vertex_id(vtxindex);
	}
      }
    }
  }
    
  // tell the primary particle copy in G4 where this output will be stored
  if (!track->GetParentID()) {
    PHG4UserPrimaryParticleInformation* userdata = static_cast<PHG4UserPrimaryParticleInformation*>(track->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
    if (userdata) {
      userdata->set_user_track_id(trackid);
      userdata->set_user_vtx_id(vtxindex);
    }
  }
  
  return;
}

void PHG4TruthTrackingAction::PostUserTrackingAction(const G4Track* track) {

  if (fpTrackingManager) {

    int trackid = track->GetTrackID();
    int primaryid = 0;
    PHG4Shower* shower = NULL;
    if ( PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()) ) {
      trackid = p->GetUserTrackId();
      primaryid = p->GetUserPrimaryId();
      shower = p->GetShower();
    }

    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries) {
      for (size_t i = 0; i < secondaries->size(); ++i) { 
	G4Track* secondary = (*secondaries)[i];    
	PHG4TrackUserInfo::SetUserParentId(const_cast<G4Track *> (secondary), trackid);
	PHG4TrackUserInfo::SetUserPrimaryId(const_cast<G4Track *> (secondary), primaryid);
	PHG4TrackUserInfo::SetShower(const_cast<G4Track *> (secondary), shower);
      }
    }
  }
  
  if ( PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()) ) {
    if ( p->GetKeep() ) {
      int trackid = p->GetUserTrackId();
      eventAction_->AddTrackidToWritelist( trackid );
    }
  }
}

void
PHG4TruthTrackingAction::SetInterfacePointers( PHCompositeNode* topNode )
{
  //now look for the map and grab a pointer to it.
  truthInfoList_ =  findNode::getClass<PHG4TruthInfoContainer>( topNode , "G4TruthInfo" );

  // if we do not find the node we need to make it.
  if ( !truthInfoList_ )
  { std::cout << "PHG4TruthEventAction::SetInterfacePointers - unable to find G4TruthInfo" << std::endl; }

}

int
PHG4TruthTrackingAction::ResetEvent(PHCompositeNode *)
{
  VertexMap.clear();
  return 0;
}
