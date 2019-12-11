#include "PHG4TruthTrackingAction.h"

#include "PHG4Particle.h"                        // for PHG4Particle
#include "PHG4Particlev2.h"
#include "PHG4Particlev3.h"
#include "PHG4Shower.h"                          // for PHG4Shower
#include "PHG4Showerv1.h"
#include "PHG4TrackUserInfoV1.h"
#include "PHG4TruthEventAction.h"
#include "PHG4TruthInfoContainer.h"
#include "PHG4UserPrimaryParticleInformation.h"
#include "PHG4VtxPointv1.h"

#include <phool/getClass.h>

#include <Geant4/G4DynamicParticle.hh>           // for G4DynamicParticle
#include <Geant4/G4ParticleDefinition.hh>        // for G4ParticleDefinition
#include <Geant4/G4PrimaryParticle.hh>
#include <Geant4/G4String.hh>                    // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Track.hh>
#include <Geant4/G4TrackingManager.hh>
#include <Geant4/G4TrackVector.hh>               // for G4TrackVector
#include <Geant4/G4VUserTrackInformation.hh>     // for G4VUserTrackInformation

#include <cmath>                                // for sqrt
#include <cstddef>                              // for size_t
#include <iostream>                              // for operator<<, endl
#include <utility>                               // for pair

using namespace std;

//________________________________________________________
PHG4TruthTrackingAction::PHG4TruthTrackingAction(PHG4TruthEventAction* eventAction)
  : m_EventAction(eventAction)
  , m_TruthInfoList(nullptr)
{
}

void PHG4TruthTrackingAction::PreUserTrackingAction(const G4Track* track)
{
  int trackid = 0;
  if (track->GetParentID())
  {
    // secondaries get negative user ids and increment downward between geant subevents
    trackid = m_TruthInfoList->mintrkindex() - 1;
  }
  else
  {
    // primaries get positive user ids and increment upward between geant subevents
    trackid = m_TruthInfoList->maxtrkindex() + 1;
  }

  // add the user id to the geant4 user info
  PHG4TrackUserInfo::SetUserTrackId(const_cast<G4Track*>(track), trackid);

  // determine the momentum vector
  G4ParticleDefinition* def = track->GetDefinition();
  int pdgid = def->GetPDGEncoding();
  double m = def->GetPDGMass();
  double ke = track->GetVertexKineticEnergy();
  double ptot = sqrt(ke * ke + 2.0 * m * ke);
  G4ThreeVector pdir = track->GetVertexMomentumDirection();
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

  ti->set_parent_id(track->GetParentID());
  if (PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()))
  {
    ti->set_parent_id(p->GetUserParentId());
  }

  ti->set_primary_id(trackid);
  if (PHG4TrackUserInfoV1* p = dynamic_cast<PHG4TrackUserInfoV1*>(track->GetUserInformation()))
  {
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

  ti->set_pid(pdgid);
  ti->set_name(def->GetParticleName());
  ti->set_e(track->GetTotalEnergy() / GeV);

  if (!track->GetParentID())
  {
    // primary track - propagate the barcode information
    PHG4UserPrimaryParticleInformation* userdata = static_cast<PHG4UserPrimaryParticleInformation*>(track->GetDynamicParticle()->GetPrimaryParticle()->GetUserInformation());
    if (userdata) ti->set_barcode(userdata->get_user_barcode());
  }

  // create a new vertex object ------------------------------------------------
  G4ThreeVector v = track->GetVertexPosition();
  map<G4ThreeVector, int>::const_iterator viter = m_VertexMap.find(v);
  int vtxindex = 0;
  if (viter != m_VertexMap.end())
  {
    vtxindex = viter->second;
  }
  else
  {
    vtxindex = m_TruthInfoList->maxvtxindex() + 1;
    if (track->GetParentID())
    {
      vtxindex = m_TruthInfoList->minvtxindex() - 1;
    }

    m_VertexMap[v] = vtxindex;
    PHG4VtxPointv1* vtxpt = new PHG4VtxPointv1(v[0] / cm,
                                               v[1] / cm,
                                               v[2] / cm,
                                               track->GetGlobalTime() / ns);
    // insert new vertex into the output
    m_TruthInfoList->AddVertex(vtxindex, vtxpt);
  }

  ti->set_vtx_id(vtxindex);

  // insert particle into the output
  m_TruthInfoList->AddParticle(trackid, ti);

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
  m_VertexMap.clear();
  return 0;
}
