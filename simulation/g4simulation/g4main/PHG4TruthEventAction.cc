#include "PHG4TruthEventAction.h"
#include "PHG4Particlev2.h"

#include "PHG4UserPrimaryParticleInformation.h"

#include "PHG4VtxPoint.h"
#include "PHG4Shower.h"

#include "PHG4TruthInfoContainer.h"
#include "PHG4HitContainer.h"
#include "PHG4Hit.h"

#include <phool/getClass.h>
#include <phool/PHPointerListIterator.h>
#include <phool/PHNode.h>

#include <Geant4/G4Event.hh>
#include <Geant4/G4TrajectoryContainer.hh>
#include <Geant4/G4VTrajectory.hh>
#include <Geant4/G4VTrajectoryPoint.hh>
#include <Geant4/G4ParticleTable.hh>
#include <Geant4/G4ParticleDefinition.hh>
#include <Geant4/globals.hh>

#include <map>

#include <Eigen/Dense>

using namespace std;

//___________________________________________________
PHG4TruthEventAction::PHG4TruthEventAction( void ):
  truthInfoList_( 0 ),
  prev_existing_lower_key( 0 ),
  prev_existing_upper_key( 0 )
{}

//___________________________________________________
void PHG4TruthEventAction::BeginOfEventAction(const G4Event* evt) {

  // if we do not find the node we need to make it.
  if ( !truthInfoList_ ) {
    std::cout << "PHG4TruthEventAction::EndOfEventAction - unable to find G4TruthInfo node" << std::endl;
    return;
  }

  const PHG4TruthInfoContainer::Map& map = truthInfoList_->GetMap();
  if (!map.empty()) {
    prev_existing_lower_key = map.begin()->first;
    prev_existing_upper_key = map.rbegin()->first;
  }
}

//___________________________________________________
void PHG4TruthEventAction::EndOfEventAction(const G4Event* evt) {
  
  // if we do not find the node we need to make it.
  if ( !truthInfoList_ ) {
    std::cout << "PHG4TruthEventAction::EndOfEventAction - unable to find G4TruthInfo node" << std::endl;
    return;
  }

  // construct a list of track ids to preserve in the the output that includes any
  // track designated in the writeList_ during processing or its ancestry chain
  
  std::set<G4int> savelist;
  std::set<int> savevtxlist;
   
  for (std::set<G4int>::const_iterator write_iter = writeList_.begin();
       write_iter != writeList_.end();
       ++write_iter) {

    std::vector<G4int> wrttracks;
    std::vector<int> wrtvtx;
    
    // usertrackid
    G4int mytrkid = *write_iter;
    PHG4Particle *particle = truthInfoList_->GetParticle(mytrkid);

    // if track is already in save list, nothing needs to be done
    if (savelist.find(mytrkid) != savelist.end()) {
      continue;
    } else {
      wrttracks.push_back(mytrkid);
      wrtvtx.push_back(particle->get_vtx_id());
    }

    // now crawl up the truth info and add parents until we hit
    // a track which is already being saved
    G4int parentid = particle->get_parent_id();
    while (savelist.find(parentid) == savelist.end() && parentid != 0) {
      particle = truthInfoList_->GetParticle(parentid);
      wrttracks.push_back(parentid);
      wrtvtx.push_back(particle->get_vtx_id());
      parentid = particle->get_parent_id();
    }
    
    // now fill the new tracks into the save list
    // while emptying the write lists
    while (wrttracks.begin() != wrttracks.end()) {
      savelist.insert(wrttracks.back());
      wrttracks.pop_back();
    }

    while (wrtvtx.begin() != wrtvtx.end()) {
      savevtxlist.insert(wrtvtx.back());
      wrtvtx.pop_back();
    }
  }
  
  // the save lists are filled now, except primary track which never
  // made it into any active volume and their vertex
  
  // loop over particles in truth list and remove them if they are not
  // in the save list and are not primary particles (parent_id == 0)
  
  int removed[4] = {0};
  PHG4TruthInfoContainer::Range truth_range = truthInfoList_->GetParticleRange();
  PHG4TruthInfoContainer::Iterator truthiter = truth_range.first;
  while (truthiter != truth_range.second) {
    
    removed[0]++;
    int trackid = truthiter->first;
    if (savelist.find(trackid) == savelist.end()) {
      // track not in save list

      // check if particle below offset 
      // primary particles had parentid = 0
      // for embedding: particles from initial sim do not have their keep flag set, we want to keep particles with trkid <= trackidoffset
      // tracks from a previous geant pass will not be recorded as leaving
      // hit in the sim, so we exclude this range from the removal
      // for regular sims, that range is zero to zero
      if (((trackid < prev_existing_lower_key)||(trackid > prev_existing_upper_key)) && ((truthiter->second)->get_parent_id() != 0)) {
	truthInfoList_->delete_particle(truthiter++);
	removed[1]++;
      } else {
	// save vertex id for primary particle which leaves no hit
	// in active area
	savevtxlist.insert((truthiter->second)->get_vtx_id());
	++truthiter;
      }
      
    } else {     
      // track was in save list, move on
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
	  truthInfoList_->AddEmbededTrkId(userdata->get_user_track_id(),userdata->get_embed());
	  truthInfoList_->AddEmbededVtxId(userdata->get_user_vtx_id(),userdata->get_embed());
	}
      }
      part = part->GetNext();
    }
    pvtx = pvtx->GetNext();
  }

  PruneShowers();
  ProcessShowers();
  
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
  
  SearchNode(topNode);
}

void PHG4TruthEventAction::SearchNode(PHCompositeNode* top) {

  // fill a lookup map between the g4hit container ids and the containers themselves
  // without knowing what the container names are in advance, only that they
  // begin G4HIT_*
  
  PHNodeIterator nodeiter(top); 
  PHPointerListIterator<PHNode> iter(nodeiter.ls());
  PHNode *thisNode; 
  while ((thisNode = iter())) {

    if (thisNode->getType() == "PHCompositeNode") {
      SearchNode(static_cast<PHCompositeNode*>(thisNode) );
    } else if (thisNode->getType() == "PHIODataNode") {
      if (thisNode->getName().find("G4HIT_") == 0) {
	PHIODataNode<PHG4HitContainer> *DNode = static_cast<PHIODataNode<PHG4HitContainer>*>(thisNode);
	if (DNode) {
	  PHG4HitContainer* object = dynamic_cast<PHG4HitContainer*>(DNode->getData());
	  if (object) {
	    hitmap_[object->GetID()] = object;
	  }
	}
      }
    }
  }
}

int PHG4TruthEventAction::ResetEvent(PHCompositeNode *) {
  writeList_.clear();
  return 0;
}

void PHG4TruthEventAction::PruneShowers() {

  PHG4TruthInfoContainer::ShowerRange range = truthInfoList_->GetShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ++iter) {
    PHG4Shower* shower = iter->second;

    std::set<int> remove_ids;
    for (PHG4Shower::ParticleIdIter jter = shower->begin_g4particle_id();
     	 jter != shower->end_g4particle_id();
     	 ++jter) {
      int g4particle_id = *jter;
      PHG4Particle* particle = truthInfoList_->GetParticle(g4particle_id);
      if (!particle) {
	remove_ids.insert(g4particle_id);
	continue;
      }
    }

    for (std::set<int>::iterator jter = remove_ids.begin();
	 jter != remove_ids.end();
	 ++jter) {
      shower->remove_g4particle_id(*jter);
    }

    std::set<int> remove_more_ids;
    for (std::map<int,std::set<PHG4HitDefs::keytype> >::iterator jter = shower->begin_g4hit_id();
	 jter != shower->end_g4hit_id();
	 ++jter) {
      int g4hitmap_id = jter->first;
      std::map<int,PHG4HitContainer*>::iterator mapiter = hitmap_.find(g4hitmap_id);
      if (mapiter == hitmap_.end()) {
	continue;
      }

      // get the g4hits from this particle in this volume
      for (std::set<PHG4HitDefs::keytype>::iterator kter = jter->second.begin();
	   kter != jter->second.end();
	   ) {
	PHG4HitDefs::keytype g4hit_id = *kter;

	PHG4Hit* g4hit = mapiter->second->findHit(g4hit_id);
	if (!g4hit) {
	  // some zero edep g4hits have been removed already
	  jter->second.erase(kter++);	  
	  continue;
	} else {
	  ++kter;
	}
      }

      if (jter->second.empty()) {
	remove_more_ids.insert(g4hitmap_id);
      }
    }

    for (std::set<int>::iterator jter = remove_more_ids.begin();
	 jter != remove_more_ids.end();
	 ++jter) {
      shower->remove_g4hit_volume(*jter);
    }    
  }

  range = truthInfoList_->GetShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ) {
    PHG4Shower* shower = iter->second;

    if (shower->empty_g4particle_id() && shower->empty_g4hit_id()) {
      if (shower ->get_edep() == 0) // check whether this shower has already been processed in the previous simulation cycles
      {
        truthInfoList_->delete_shower(iter++);
        continue;
      }
    }

    ++iter;
  }
  
}

void PHG4TruthEventAction::ProcessShowers() {

  PHG4TruthInfoContainer::ShowerRange range = truthInfoList_->GetShowerRange();
  for (PHG4TruthInfoContainer::ShowerIterator iter = range.first;
       iter != range.second;
       ++iter) {
    PHG4Shower* shower = iter->second;
    
    // Data structures to hold weighted pca
    std::vector<std::vector<float> > points;
    std::vector<float> weights;
    float sumw = 0.0;
    float sumw2 = 0.0;

    for (std::map<int,std::set<PHG4HitDefs::keytype> >::iterator iter = shower->begin_g4hit_id();
	 iter != shower->end_g4hit_id();
	 ++iter) {
      int g4hitmap_id = iter->first;
      std::map<int,PHG4HitContainer*>::iterator mapiter = hitmap_.find(g4hitmap_id);
      if (mapiter == hitmap_.end()) {
	continue;
      }

      PHG4HitContainer* hits = mapiter->second;

      unsigned int nhits = 0;
      float edep = 0.0;
      float eion = 0.0;
      float light_yield = 0.0;
      float edep_e = 0.0;
      float edep_h = 0.0;

      // get the g4hits from this particle in this volume
      for (std::set<PHG4HitDefs::keytype>::iterator kter = iter->second.begin();
	   kter != iter->second.end();
	   ++kter) {
	PHG4HitDefs::keytype g4hit_id = *kter;

	PHG4Hit* g4hit = hits->findHit(g4hit_id);
	if (!g4hit) {
	  cout << PHWHERE << " missing g4hit" << endl;
	  continue;
	}

	PHG4Particle* particle = truthInfoList_->GetParticle(g4hit->get_trkid());
	if (!particle) {	  
	  cout << PHWHERE << " missing g4particle for track " 
	       << g4hit->get_trkid() << endl;
	  continue;
	}

	PHG4VtxPoint* vtx = truthInfoList_->GetVtx( particle->get_vtx_id() );
	if (!vtx) {
	  cout << PHWHERE << " missing g4vertex" << endl;
	  continue;
	}
	
	// shower location and shape info
	
	if (!isnan(g4hit->get_x(0)) &&
	    !isnan(g4hit->get_y(0)) &&
	    !isnan(g4hit->get_z(0))) {
	  
	  std::vector<float> entry(3);
	  entry[0] = g4hit->get_x(0);
	  entry[1] = g4hit->get_y(0);
	  entry[2] = g4hit->get_z(0);
	  
	  points.push_back(entry);
	  float w = g4hit->get_edep();
	  weights.push_back(w);
	  sumw += w;
	  sumw2 += w*w;
	}

	if (!isnan(g4hit->get_x(1)) &&
	    !isnan(g4hit->get_y(1)) &&
	    !isnan(g4hit->get_z(1))) {
	  
	  std::vector<float> entry(3);
	  entry[0] = g4hit->get_x(1);
	  entry[1] = g4hit->get_y(1);
	  entry[2] = g4hit->get_z(1);
	  
	  points.push_back(entry);	  
	  float w = g4hit->get_edep();
	  weights.push_back(w);
	  sumw += w;
	  sumw2 += w*w;
	}
	
	// e/h ratio
	
	if (!isnan(g4hit->get_edep())) {
	  if (abs(particle->get_pid()) == 11) {
	    edep_e += g4hit->get_edep();
	  } else {
	    edep_h += g4hit->get_edep();
	  }
	}

	// summary info
	
	if (g4hit)                                   ++nhits;
	if (!isnan(g4hit->get_edep()))               edep += g4hit->get_edep();
	if (!isnan(g4hit->get_eion()))               eion += g4hit->get_eion();
	if (!isnan(g4hit->get_light_yield())) light_yield += g4hit->get_light_yield();
      } // g4hit loop

      // summary info
      
      if (nhits)              shower->set_nhits(g4hitmap_id,nhits);
      if (edep != 0.0)        shower->set_edep(g4hitmap_id,edep);
      if (eion != 0.0)        shower->set_eion(g4hitmap_id,eion);
      if (light_yield != 0.0) shower->set_light_yield(g4hitmap_id,light_yield);
      if (edep_h != 0.0)      shower->set_eh_ratio(g4hitmap_id,edep_e/edep_h);
    } // volume loop

    // fill Eigen matrices to compute wPCA
    // resizing these non-destructively is expensive
    // so I fill vectors and then copy
    Eigen::Matrix<double, Eigen::Dynamic, 3> X(points.size(),3);
    Eigen::Matrix<double, Eigen::Dynamic, 1> W(weights.size(),1);

    for (unsigned int i=0; i<points.size(); ++i) {
      for (unsigned int j=0; j<3; ++j)  {
	X(i,j) = points[i][j];
      }
      W(i,0) = weights[i];
    }

    // mean value of shower
    double prefactor = 1.0 / sumw;
    Eigen::Matrix<double, 1, 3> mean = prefactor * W.transpose() * X;

    // compute residual relative to the mean
    for (unsigned int i=0; i<points.size(); ++i) {
      for (unsigned int j=0; j<3; ++j) X(i,j) = points[i][j] - mean(0,j);
    }

    // weighted covariance matrix
    prefactor = sumw / (pow(sumw,2) - sumw2); // effectivelly 1/(N-1) when w_i = 1.0
    Eigen::Matrix<double, 3, 3> covar = prefactor * (X.transpose() * W.asDiagonal() * X);
       
    shower->set_x(mean(0,0));
    shower->set_y(mean(0,1));
    shower->set_z(mean(0,2));

    for (unsigned int i = 0; i < 3; ++i) {
      for (unsigned int j = 0; j <= i; ++j) {
	shower->set_covar(i,j,covar(i,j));
      }
    }

    // shower->identify();
  }
}
