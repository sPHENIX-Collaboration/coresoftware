#include "PHG4TruthTrackingAction.h"

#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Shower.h"  // for PHG4Shower
#include "PHG4Showerv1.h"
#include "PHG4TrackUserInfoV1.h"
#include "PHG4TruthEventAction.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4UserPrimaryParticleInformation.h"

#include <phool/getClass.h>

#include <Geant4/G4PrimaryParticle.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackingManager.hh>
#include <Geant4/G4VUserTrackInformation.hh>  // for G4VUserTrackInformation

#include <cmath>     // for sqrt
#include <cstddef>   // for size_t
#include <iostream>  // for operator<<, endl
#include <utility>   // for pair

using namespace std;

//________________________________________________________
PHG4TruthTrackingAction::PHG4TruthTrackingAction(PHG4TruthEventAction* eventAction)
  : m_EventAction(eventAction)
  , m_TruthInfoList(nullptr)
{
}

void PHG4TruthTrackingAction::PreUserTrackingAction(const G4Track* track)
{
  // insert particle into the output
  PHG4Particle* ti = (*m_TruthInfoList->AddParticle(const_cast<G4Track*>(track))).second;
  int trackid = ti->get_track_id();
  int vtxindex = ti->get_vtx_id();

  // create or add to a new shower object --------------------------------------
  if (!track->GetParentID())
  {
    PHG4Showerv1* shower = new PHG4Showerv1();
    PHG4TrackUserInfo::SetShower(const_cast<G4Track*>(track), shower);
    m_TruthInfoList->AddShower(trackid, shower);
    shower->set_id(trackid);  // fyi, secondary showers may not share these ids
    shower->set_parent_particle_id(trackid);
    shower->set_parent_shower_id(0);
  }
  else
  {
    // get shower
    if (G4VUserTrackInformation* p = track->GetUserInformation())
    {
      if (PHG4TrackUserInfoV1* pp = dynamic_cast<PHG4TrackUserInfoV1*>(p))
      {
        if (pp->GetShower())
        {
          pp->GetShower()->add_g4particle_id(trackid);
          pp->GetShower()->add_g4vertex_id(vtxindex);
        }
      }
    }
  }

  // tell the primary particle copy in G4 where this output will be stored
  if (!track->GetParentID())
  {
    PHG4UserPrimaryParticleInformation* userdata = static_cast<PHG4UserPrimaryParticleInformation*>(track->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
    if (userdata)
    {
      userdata->set_user_track_id(trackid);
      userdata->set_user_vtx_id(vtxindex);
    }
  }

  return;
}

void PHG4TruthTrackingAction::PostUserTrackingAction(const G4Track* track)
{
  if (fpTrackingManager)
  {
    int trackid = track->GetTrackID();
    int primaryid = 0;
    PHG4Shower* shower = nullptr;
    if (PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()))
    {
      trackid = p->GetUserTrackId();
      primaryid = p->GetUserPrimaryId();
      shower = p->GetShower();
    }

    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries)
    {
      for (size_t i = 0; i < secondaries->size(); ++i)
      {
        G4Track* secondary = (*secondaries)[i];
        PHG4TrackUserInfo::SetUserParentId(const_cast<G4Track*>(secondary), trackid);
        PHG4TrackUserInfo::SetUserPrimaryId(const_cast<G4Track*>(secondary), primaryid);
        PHG4TrackUserInfo::SetShower(const_cast<G4Track*>(secondary), shower);
      }
    }
  }

  if (PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()))
  {
    if (p->GetKeep())
    {
      int trackid = p->GetUserTrackId();
      m_EventAction->AddTrackidToWritelist(trackid);
    }
  }
}

void PHG4TruthTrackingAction::SetInterfacePointers(PHCompositeNode* topNode)
{
  //now look for the map and grab a pointer to it.
  m_TruthInfoList = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // if we do not find the node we need to make it.
  if (!m_TruthInfoList)
  {
    std::cout << "PHG4TruthEventAction::SetInterfacePointers - unable to find G4TruthInfo" << std::endl;
  }
}

int PHG4TruthTrackingAction::ResetEvent(PHCompositeNode*)
{
  return 0;
}
