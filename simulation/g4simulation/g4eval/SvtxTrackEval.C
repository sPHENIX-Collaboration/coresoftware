
#include "SvtxTrackEval.h"

#include "SvtxClusterEval.h"

#include <fun4all/getClass.h>
#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>

#include <cstdlib>
#include <set>
#include <float.h>
#include <algorithm>

using namespace std;

SvtxTrackEval::SvtxTrackEval(PHCompositeNode* topNode)
  : _topNode(topNode) {
}

std::set<PHG4Hit*> SvtxTrackEval::all_truth_hits(SvtxTrack* track) {

  // need things off of the DST...
  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(_topNode,"SvtxClusterMap");
  if (!clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }
  
  std::set<PHG4Hit*> truth_hits;

  // loop over all clusters...
  for (unsigned int ilayer = 0; ilayer < 100; ++ilayer) {
    if (!track->hasCluster(ilayer)) continue;

    SvtxCluster* cluster = clustermap->get(track->getClusterID(ilayer));
    SvtxClusterEval clustereval(_topNode);

    std::set<PHG4Hit*> new_hits = clustereval.all_truth_hits(cluster);

    std::set<PHG4Hit*> union_hits; // placeholder for union of new hits and truth hits

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   new_hits.begin(),new_hits.end(),
		   std::inserter(union_hits,union_hits.begin()));

    std::swap(truth_hits,union_hits); // swap union into truth_hits
  }

  return truth_hits;
}
  
std::set<PHG4Particle*> SvtxTrackEval::all_truth_particles(SvtxTrack* track) {

  // need things off of the DST...
  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(_topNode,"SvtxClusterMap");
  if (!clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }
  
  std::set<PHG4Particle*> truth_particles;

  // loop over all clusters...
  for (unsigned int ilayer = 0; ilayer < 100; ++ilayer) {
    if (!track->hasCluster(ilayer)) continue;

    SvtxCluster* cluster = clustermap->get(track->getClusterID(ilayer));
    SvtxClusterEval clustereval(_topNode);

    std::set<PHG4Particle*> new_particles = clustereval.all_truth_particles(cluster);

    std::set<PHG4Particle*> union_particles; // placeholder for union of new particles and truth particles

    std::set_union(truth_particles.begin(),truth_particles.end(),
		   new_particles.begin(),new_particles.end(),
		   std::inserter(union_particles,union_particles.begin()));

    std::swap(truth_particles,union_particles); // swap union into truth_particles
  }

  return truth_particles;
}

PHG4Particle* SvtxTrackEval::max_truth_particle_by_nclusters(SvtxTrack* track) {

  std::set<PHG4Particle*> particles = all_truth_particles(track);

  PHG4Particle* max_particle = NULL;
  unsigned int max_nclusters = 0;
  
  for(std::set<PHG4Particle*>::iterator iter = particles.begin();
      iter != particles.end();
      ++iter) {
    PHG4Particle* candidate = *iter;
    unsigned int nclusters = get_nclusters_contribution(track,candidate);
    if (max_nclusters > nclusters) {
      max_nclusters = nclusters;
      max_particle = candidate;
    }
  }
  
  return max_particle;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Particle* truthparticle) { 

  // need things off of the DST...
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(_topNode,"SvtxTrackMap");
  if (!trackmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
    exit(-1);
  }

  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(_topNode,"SvtxClusterMap");
  if (!clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }

  std::set<SvtxTrack*> tracks;
  
  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter) {
    SvtxTrack* track = &iter->second;
    
    // loop over all clusters    
    for (unsigned int ilayer = 0; ilayer < 100; ++ilayer) {
      if (!track->hasCluster(ilayer)) continue;

      SvtxCluster* cluster = clustermap->get(track->getClusterID(ilayer));
      SvtxClusterEval clustereval(_topNode);

      // loop over all particles
      std::set<PHG4Particle*> particles = clustereval.all_truth_particles(cluster);
      for (std::set<PHG4Particle*>::iterator jter = particles.begin();
	   jter != particles.end();
	   ++jter) {
	PHG4Particle* candidate = *jter;
	// if track id matches argument add to output
	if (candidate->get_track_id() == truthparticle->get_track_id()) {
	  tracks.insert(track);
	}
      }
    }
  }

  return tracks;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Hit* truthhit) {

  // need things off of the DST...
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(_topNode,"SvtxTrackMap");
  if (!trackmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
    exit(-1);
  }

  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(_topNode,"SvtxClusterMap");
  if (!clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }

  std::set<SvtxTrack*> tracks;
  
  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter) {
    SvtxTrack* track = &iter->second;
    
    // loop over all clusters    
    for (unsigned int ilayer = 0; ilayer < 100; ++ilayer) {
      if (!track->hasCluster(ilayer)) continue;

      SvtxCluster* cluster = clustermap->get(track->getClusterID(ilayer));
      SvtxClusterEval clustereval(_topNode);

      // loop over all hits
      std::set<PHG4Hit*> hits = clustereval.all_truth_hits(cluster);
      for (std::set<PHG4Hit*>::iterator jter = hits.begin();
	   jter != hits.end();
	   ++jter) {
	PHG4Hit* candidate = *jter;
	// if track id matches argument add to output
	if (candidate->get_trkid() == truthhit->get_trkid()) {
	  tracks.insert(track);
	}
      }
    }
  }

  return tracks;
}
  
// overlap calculations
unsigned int SvtxTrackEval::get_nclusters_contribution(SvtxTrack* track, PHG4Particle* particle) {

  SvtxClusterMap* clustermap = findNode::getClass<SvtxClusterMap>(_topNode,"SvtxClusterMap");
  if (!clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }
  
  unsigned int nclusters = 0;

  // loop over all clusters    
  for (unsigned int ilayer = 0; ilayer < 100; ++ilayer) {
    if (!track->hasCluster(ilayer)) continue;

    SvtxCluster* cluster = clustermap->get(track->getClusterID(ilayer));
    SvtxClusterEval clustereval(_topNode);

    // loop over all particles
    std::set<PHG4Particle*> particles = clustereval.all_truth_particles(cluster);
    for (std::set<PHG4Particle*>::iterator jter = particles.begin();
	 jter != particles.end();
	 ++jter) {
      PHG4Particle* candidate = *jter;
      // if track id matches argument add to output
      if (candidate->get_track_id() == particle->get_track_id()) {
	++nclusters;
      }
    }
  }
  
  return nclusters;
}
