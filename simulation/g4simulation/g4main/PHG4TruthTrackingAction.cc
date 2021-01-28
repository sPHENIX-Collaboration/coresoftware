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
  , m_G4ParticleStack()
{
}

void PHG4TruthTrackingAction::PreUserTrackingAction(const G4Track* track)
{
  // insert particle into the output
  PHG4Particle* ti = (*m_TruthInfoList->AddParticle(const_cast<G4Track*>(track))).second;

  // Negative G4 track id values indicate unwanted tracks to be deleted
  // Initially all tracks except primary ones flagged as unwanted
  int track_id_g4 = track->GetTrackID() * (track->GetParentID() ? -1 : +1);
  int trackid = ti->get_track_id();
  int vtxindex = ti->get_vtx_id();

  m_CurrG4Particle = {track_id_g4, trackid, vtxindex};

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

  UpdateG4ParticleStack(track);
}

/**
 * Updates the stack of parent particles and removes unwanted ones from the
 * truth info container.
 */
void PHG4TruthTrackingAction::UpdateG4ParticleStack(const G4Track* track)
{
  while (!m_G4ParticleStack.empty())
  {
    if ( std::abs(m_G4ParticleStack.back().g4track_id) == track->GetParentID() )
    {
      break;
    }
    else
    {
      if (m_G4ParticleStack.back().g4track_id < 0)
      {
        m_TruthInfoList->delete_particle( m_G4ParticleStack.back().particle_id );
      }
      m_G4ParticleStack.pop_back();
    }
  }

  m_G4ParticleStack.push_back(m_CurrG4Particle);

  // Change sign of G4 track id of all upstream tracks in the stack to positive
  // in order to keep the track
  PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation());
  bool keep_curr_track = p && p->GetKeep() ? true : false;

  auto stack_iter = m_G4ParticleStack.rbegin();
  while (keep_curr_track && stack_iter != m_G4ParticleStack.rend() && stack_iter->g4track_id < 0)
  {
    stack_iter->g4track_id *= -1;
    ++stack_iter;
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
  while (!m_G4ParticleStack.empty())
  {
    if ( m_G4ParticleStack.back().g4track_id < 0)
    {
      m_TruthInfoList->delete_particle( m_G4ParticleStack.back().particle_id );
    }
    m_G4ParticleStack.pop_back();
  }

  return 0;
}
