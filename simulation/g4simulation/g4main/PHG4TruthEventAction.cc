#include "PHG4TruthEventAction.h"
#include "PHG4Particlev2.h"

#include "PHG4UserPrimaryParticleInformation.h"

#include "PHG4VtxPoint.h"
#include "PHG4TruthInfoContainer.h"

#include <phool/getClass.h>

#include <Geant4/G4Event.hh>
#include <Geant4/G4TrajectoryContainer.hh>
#include <Geant4/G4VTrajectory.hh>
#include <Geant4/G4VTrajectoryPoint.hh>
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/globals.hh>


#include <map>

using namespace std;

//___________________________________________________
PHG4TruthEventAction::PHG4TruthEventAction( void ):
  truthInfoList_( 0 )
{}

//___________________________________________________
void PHG4TruthEventAction::BeginOfEventAction(const G4Event* evt) {}

//___________________________________________________
void PHG4TruthEventAction::EndOfEventAction(const G4Event* evt)
{
  // if we do not find the node we need to make it.
  if ( !truthInfoList_ )
    {
      std::cout << "PHG4TruthEventAction::EndOfEventAction - unable to find G4TruthInfo node" << std::endl;
      return;
    }
  
  set<G4int> savelist;
  set<int> savevtxlist;
  set<G4int>::const_iterator write_iter;
  vector<G4int> wrttracks;
  vector<int> wrtvtx;
  for (write_iter = writeList_.begin(); write_iter != writeList_.end(); ++write_iter)
    {
      G4int mytrkid = *write_iter;
      PHG4Particle *particle = truthInfoList_->GetParticle(mytrkid);
      // if track is already in save list, nothing needs to be done
      if (savelist.find(mytrkid) != savelist.end())
        {
          continue;
        }
      else
        {
          wrttracks.push_back(mytrkid);
          wrtvtx.push_back(particle->get_vtx_id());
        }
      // now crawl up the truth info and add parents until we hit
      // a track which is already being saved
      G4int parentid = particle->get_parent_id();
      while (savelist.find(parentid) == savelist.end() && parentid > 0)
        {
          particle = truthInfoList_->GetParticle(parentid);
          wrttracks.push_back(parentid);
          wrtvtx.push_back(particle->get_vtx_id());
          parentid = particle->get_parent_id();
        }
      // now fill the new tracks into the save list
      while (wrttracks.begin() != wrttracks.end())
        {
          savelist.insert(wrttracks.back());
          wrttracks.pop_back();
        }
      while (wrtvtx.begin() != wrtvtx.end())
        {
          savevtxlist.insert(wrtvtx.back());
          wrtvtx.pop_back();
        }
    }
  
  // the save lists are filled now, except primary track which never
  // made it into any active volume and their vertex
  // loop over particles in truth list and remove them if they are not
  // in the save list andare not primary particles (parent_id == 0)
  int removed[4] = {0};
  PHG4TruthInfoContainer::Range truth_range = truthInfoList_->GetParticleRange();
  PHG4TruthInfoContainer::Iterator truthiter = truth_range.first;
  while (truthiter != truth_range.second)
    {
      removed[0]++;
      int trackid = truthiter->first;
      if (savelist.find(trackid) == savelist.end()) // track not in save list
        {
          // check if particle below offset 
	  // primary particles had parentid = 0
	  // for embedding: particles from initial sim do not have their keep flag set, we want to keep particles with trkid <= trackidoffset
	  // trackidoffset is the geantid of the last particle
	  // for regular sims, trackidoffset is zero
          if ((truthiter->second)->get_parent_id() != 0 && truthiter->first > 0) // MPM I'll need to revise this section
            {
              truthInfoList_->delete_particle(truthiter++);
              removed[1]++;
            }
          else
            {
              // save vertex id for primary particle which leaves no hit
              // in active area
              savevtxlist.insert((truthiter->second)->get_vtx_id());
              ++truthiter;
            }
        }
      else
        {
          ++truthiter;
        }
    }
  PHG4TruthInfoContainer::VtxRange vtxrange = truthInfoList_->GetVtxRange();
  PHG4TruthInfoContainer::VtxIterator vtxiter = vtxrange.first;
  while (vtxiter != vtxrange.second)
    {
      removed[2]++;
      if (savevtxlist.find(vtxiter->first) == savevtxlist.end())
        {
          truthInfoList_->delete_vtx(vtxiter++);
          removed[3]++;
        }
      else
        {
          ++vtxiter;
        }
    }

  // loop over all input particles and fish out the ones which have the embed flag set
  // and store their geant track ids in truthinfo container
  G4PrimaryVertex *pvtx = evt->GetPrimaryVertex();
  while (pvtx) {

    G4PrimaryParticle *part = pvtx->GetPrimary();
    while (part) {

      PHG4UserPrimaryParticleInformation *userdata = dynamic_cast<PHG4UserPrimaryParticleInformation *> (part->GetUserInformation());
      if (userdata) {
	if (userdata->get_embed()) {
	  truthInfoList_->AddEmbededTrkId(userdata->get_user_track_id()); // use primary ID for the embedded list
	}
      }
      part = part->GetNext();
    }
    pvtx = pvtx->GetNext();
  }
  
  return;
}

//___________________________________________________
void PHG4TruthEventAction::AddTrackidToWritelist(const G4int trackid) {
   writeList_.insert(trackid);
}

//___________________________________________________
void PHG4TruthEventAction::SetInterfacePointers(PHCompositeNode* topNode) {
  
  //now look for the map and grab a pointer to it.
  truthInfoList_ =  findNode::getClass<PHG4TruthInfoContainer>( topNode , "G4TruthInfo" );

  // if we do not find the node we need to make it.
  if ( !truthInfoList_ ) {
    std::cout << "PHG4TruthEventAction::SetInterfacePointers - unable to find G4TruthInfo" << std::endl;
  }
}

int PHG4TruthEventAction::ResetEvent(PHCompositeNode *) {
  writeList_.clear();
  return 0;
}
