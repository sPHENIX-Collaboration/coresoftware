#include "PHG4TruthTrackingAction.h"
#include "PHG4TruthEventAction.h"
#include <PHG4TruthInfoContainer.h>
#include "PHG4TrackUserInfoV1.h"
#include "PHG4Particlev2.h"
#include "PHG4VtxPointv1.h"

#include <phool/getClass.h>

#include <Geant4/G4Step.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackingManager.hh>
#include <Geant4/G4TrackVector.hh>

using namespace std;

const int VERBOSE = 0;

//________________________________________________________
PHG4TruthTrackingAction::PHG4TruthTrackingAction( PHG4TruthEventAction* eventAction ) :
  primarytrackidoffset(0),
  secondarytrackidoffset(0),
  eventAction_( eventAction ), 
  truthInfoList_( NULL )
{}

void
PHG4TruthTrackingAction::PreUserTrackingAction( const G4Track* track)
{
  G4ThreeVector v = track->GetVertexPosition();
  G4ThreeVector pdir = track->GetVertexMomentumDirection();
  int trackid = 0;
  // if (track->GetParentID())
  //   {
  //     trackid = -track->GetTrackID() + secondarytrackidoffset;
  //   }
  // else
  //   {
       trackid = track->GetTrackID() + primarytrackidoffset;
  //   }
  PHG4TrackUserInfo::SetUserTrackId(const_cast<G4Track *> (track), trackid);

  
  G4ParticleDefinition* def = track->GetDefinition();
  int pdgid = def->GetPDGEncoding();
  //   double charge = def->GetPDGCharge();
  double m = def->GetPDGMass();
  double ke = track->GetVertexKineticEnergy();
  double ptot = sqrt(ke * ke + 2.0 * m * ke);
  pdir *= ptot;

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
  map<G4ThreeVector, int>::const_iterator viter = VertexMap.find(v);
  int vtxindex;
  if (viter == VertexMap.end())
    {
      vtxindex = truthInfoList_->maxvtxindex() + 1;
      VertexMap[v] = vtxindex;
      PHG4VtxPointv1 *vtxpt = new PHG4VtxPointv1(v[0] / cm, v[1] / cm, v[2] / cm, track->GetGlobalTime() / ns);
      // 	  cout << "adding vertex " << vtxindex << endl;
      //  	  vtxpt->identify();
      truthInfoList_->AddVertex(vtxindex, vtxpt);
    }
  else
    {
      vtxindex = viter->second;
      // 	  cout << "found vertex " << vtxindex << endl;
    }
  ti->set_vtx_id(vtxindex);
  //       cout << "Adding particle trkid: " << trackid << endl;
  //       ti->identify();
  truthInfoList_->AddParticle(trackid, ti);
  return;
}

void PHG4TruthTrackingAction::PostUserTrackingAction(const G4Track* track) {

  if (fpTrackingManager) {

    int trackid = track->GetTrackID();
    int primaryid = 0;
    if ( PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()) ) {
      trackid = p->GetUserTrackId();
      primaryid = p->GetUserPrimaryId();
    }

    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries) {
      for (size_t i = 0; i < secondaries->size(); ++i) { 
	G4Track* secondary = (*secondaries)[i];    
	PHG4TrackUserInfo::SetUserParentId(const_cast<G4Track *> (secondary), trackid);
	PHG4TrackUserInfo::SetUserPrimaryId(const_cast<G4Track *> (secondary), primaryid);
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
