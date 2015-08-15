
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
  : _topNode(topNode),
    _clustereval(topNode),
    _cache_all_truth_hits(),
    _cache_all_truth_particles(),
    _cache_max_truth_particle_by_nclusters(),
    _cache_all_tracks_from_particle(),
    _cache_best_track_from_particle(),
    _cache_all_tracks_from_g4hit(),
    _cache_get_nclusters_contribution() {
}


void SvtxTrackEval::next_event(PHCompositeNode* topNode) {

  _cache_all_truth_hits.clear();
  _cache_all_truth_particles.clear();
  _cache_max_truth_particle_by_nclusters.clear();
  _cache_all_tracks_from_particle.clear();
  _cache_best_track_from_particle.clear();
  _cache_all_tracks_from_g4hit.clear();
  _cache_get_nclusters_contribution.clear();
  
  _clustereval.next_event(topNode);
  
  _topNode = topNode;  
}

std::set<PHG4Hit*> SvtxTrackEval::all_truth_hits(SvtxTrack* track) {

  if (_cache_all_truth_hits.find(track) != _cache_all_truth_hits.end()) {
    return _cache_all_truth_hits[track];
  }
  
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

    std::set<PHG4Hit*> new_hits = _clustereval.all_truth_hits(cluster);

    std::set<PHG4Hit*> union_hits; // placeholder for union of new hits and truth hits

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   new_hits.begin(),new_hits.end(),
		   std::inserter(union_hits,union_hits.begin()));

    std::swap(truth_hits,union_hits); // swap union into truth_hits
  }

  _cache_all_truth_hits.insert(make_pair(track,truth_hits));
  
  return truth_hits;
}
  
std::set<PHG4Particle*> SvtxTrackEval::all_truth_particles(SvtxTrack* track) {

  if (_cache_all_truth_particles.find(track) != _cache_all_truth_particles.end()) {
    return _cache_all_truth_particles[track];
  }
  
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

    std::set<PHG4Particle*> new_particles = _clustereval.all_truth_particles(cluster);

    std::set<PHG4Particle*> union_particles; // placeholder for union of new particles and truth particles

    std::set_union(truth_particles.begin(),truth_particles.end(),
		   new_particles.begin(),new_particles.end(),
		   std::inserter(union_particles,union_particles.begin()));

    std::swap(truth_particles,union_particles); // swap union into truth_particles
  }

  _cache_all_truth_particles.insert(make_pair(track,truth_particles));

  return truth_particles;
}

PHG4Particle* SvtxTrackEval::max_truth_particle_by_nclusters(SvtxTrack* track) {

  if (_cache_max_truth_particle_by_nclusters.find(track) !=
      _cache_max_truth_particle_by_nclusters.end()) {
    return _cache_max_truth_particle_by_nclusters[track];
  }
  
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

  _cache_max_truth_particle_by_nclusters.insert(make_pair(track,max_particle));
  
  return max_particle;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Particle* truthparticle) { 

  if (_cache_all_tracks_from_particle.find(truthparticle) !=
      _cache_all_tracks_from_particle.end()) {
    return _cache_all_tracks_from_particle[truthparticle];
  }
  
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

      // loop over all particles
      std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster);
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

  _cache_all_tracks_from_particle.insert(make_pair(truthparticle,tracks));
  
  return tracks;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Hit* truthhit) {

  if (_cache_all_tracks_from_g4hit.find(truthhit) !=
      _cache_all_tracks_from_g4hit.end()) {
    return _cache_all_tracks_from_g4hit[truthhit];
  }
  
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

      // loop over all hits
      std::set<PHG4Hit*> hits = _clustereval.all_truth_hits(cluster);
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

  _cache_all_tracks_from_g4hit.insert(make_pair(truthhit,tracks));

  return tracks;
}

SvtxTrack* SvtxTrackEval::best_track_from(PHG4Particle* truthparticle) { 

  if (_cache_best_track_from_particle.find(truthparticle) !=
      _cache_best_track_from_particle.end()) {
    return _cache_best_track_from_particle[truthparticle];
  }

  SvtxTrack* best_track = NULL;
  unsigned int best_purity = 0;
  std::set<SvtxTrack*> tracks = all_tracks_from(truthparticle);
  for (std::set<SvtxTrack*>::iterator iter = tracks.begin();
       iter != tracks.end();
       ++iter) {
    SvtxTrack* track = *iter;
    unsigned int purity = get_nclusters_contribution(track,truthparticle);
    if (purity > best_purity) {
      best_track = track;
      best_purity = purity;
    }
  }
  
  _cache_best_track_from_particle.insert(make_pair(truthparticle,best_track));
  
  return best_track;
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

    // loop over all particles
    std::set<PHG4Particle*> particles = _clustereval.all_truth_particles(cluster);
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
  
  _cache_get_nclusters_contribution.insert(make_pair(make_pair(track,particle),nclusters));
  
  return nclusters;
}
