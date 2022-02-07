#include "PHG4TruthTrackingAction.h"

#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Particlev2.h"
#include "PHG4Particlev3.h"
#include "PHG4Shower.h"  // for PHG4Shower
#include "PHG4Showerv1.h"
#include "PHG4TrackUserInfoV1.h"
#include "PHG4TruthEventAction.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4UserPrimaryParticleInformation.h"
#include "PHG4VtxPointv1.h"

#include <phool/getClass.h>

#include <Geant4/G4DynamicParticle.hh>     // for G4DynamicParticle
#include <Geant4/G4ParticleDefinition.hh>  // for G4ParticleDefinition
#include <Geant4/G4PrimaryParticle.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackVector.hh>  // for G4TrackVector
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
  , m_CurrG4Particle()
{
}

void PHG4TruthTrackingAction::PreUserTrackingAction(const G4Track* track)
{
  // insert particle into the output
  PHG4Particle* ti = AddParticle(*m_TruthInfoList, *const_cast<G4Track*>(track));

  // Negative G4 track id values indicate unwanted tracks to be deleted
  // Initially all tracks except primary ones flagged as unwanted
  int track_id_g4 = track->GetTrackID() * (track->GetParentID() ? -1 : +1);
  int trackid = ti->get_track_id();

  // add the user id to the geant4 user info
  PHG4TrackUserInfo::SetUserTrackId(const_cast<G4Track*>(track), trackid);

  if (PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()))
  {
    ti->set_parent_id(p->GetUserParentId());

    if (track->GetParentID())
    {
      ti->set_primary_id(p->GetUserPrimaryId());
    }
    else
    {
      PHG4TrackUserInfo::SetUserPrimaryId(const_cast<G4Track*>(track), trackid);
      ti->set_primary_id(trackid);
    }
  }

  if (!track->GetParentID())
  {
    // primary track - propagate the barcode information
    PHG4UserPrimaryParticleInformation* userdata = static_cast<PHG4UserPrimaryParticleInformation*>(track->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
    if (userdata) ti->set_barcode(userdata->get_user_barcode());
  }

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
  m_VertexMap.clear();

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

PHG4Particle* PHG4TruthTrackingAction::AddParticle(PHG4TruthInfoContainer& truth, G4Track& track)
{
  int trackid = 0;
  if (track.GetParentID())
  {
    // secondaries get negative user ids and increment downward between geant subevents
    trackid = truth.mintrkindex() - 1;
  }
  else
  {
    // primaries get positive user ids and increment upward between geant subevents
    trackid = truth.maxtrkindex() + 1;
  }

  // determine the momentum vector
  G4ParticleDefinition* def = track.GetDefinition();
  int pdgid = def->GetPDGEncoding();
  double mass = def->GetPDGMass();
  double ke = track.GetVertexKineticEnergy();
  double ptot = sqrt(ke * ke + 2.0 * mass * ke);
  G4ThreeVector pdir = track.GetVertexMomentumDirection();
  pdir *= ptot;
  PHG4Particle* ti = nullptr;
  // create a new particle -----------------------------------------------------
  if (def->IsGeneralIon())  // for ions save a and z in v3 of phg4particle
  {
    ti = new PHG4Particlev3();
    ti->set_A(def->GetAtomicMass());
    ti->set_Z(def->GetAtomicNumber());
  }
  else
  {
    ti = new PHG4Particlev2;
  }
  ti->set_px(pdir[0] / GeV);
  ti->set_py(pdir[1] / GeV);
  ti->set_pz(pdir[2] / GeV);
  ti->set_track_id(trackid);

  ti->set_parent_id(track.GetParentID());
  ti->set_primary_id(trackid);

  ti->set_pid(pdgid);
  ti->set_name(def->GetParticleName());
  ti->set_e(track.GetTotalEnergy() / GeV);

  // Add new or reuse a vertex and let the track know about it
  PHG4VtxPoint* vtx = AddVertex(truth, track);
  ti->set_vtx_id(vtx->get_id());

  return truth.AddParticle(trackid, ti)->second;
}

PHG4VtxPoint* PHG4TruthTrackingAction::AddVertex(PHG4TruthInfoContainer& truth, const G4Track& track)
{
  G4ThreeVector v = track.GetVertexPosition();
  int vtxindex = (track.GetParentID() == 0 ? truth.maxvtxindex() + 1 : truth.minvtxindex() - 1);

  auto [iter, inserted] = m_VertexMap.insert(std::make_pair(v, vtxindex));

  // If could not add a unique vertex => return the existing one
  if (!inserted)
  {
    return truth.GetVtxMap().find(iter->second)->second;
  }
  // otherwise, create and add a new one
  PHG4VtxPoint* vtxpt = new PHG4VtxPointv1(v[0]/cm, v[1]/cm, v[2]/cm, track.GetGlobalTime()/ns, vtxindex);

  return truth.AddVertex(vtxindex, vtxpt)->second;
}
