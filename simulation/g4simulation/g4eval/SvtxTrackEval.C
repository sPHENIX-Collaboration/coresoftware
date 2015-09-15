
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
  : _clustereval(topNode),
    _trackmap(NULL),
    _clustermap(NULL),
    _do_cache(true),
    _cache_all_truth_hits(),
    _cache_all_truth_particles(),
    _cache_max_truth_particle_by_nclusters(),
    _cache_all_tracks_from_particle(),
    _cache_best_track_from_particle(),
    _cache_all_tracks_from_g4hit(),
    _cache_get_nclusters_contribution() {
  get_node_pointers(topNode);
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
  
  get_node_pointers(topNode);
}

std::set<PHG4Hit*> SvtxTrackEval::all_truth_hits(SvtxTrack* track) {

  if (_do_cache) {
    std::map<SvtxTrack*,std::set<PHG4Hit*> >::iterator iter =
      _cache_all_truth_hits.find(track);
    if (iter != _cache_all_truth_hits.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Hit*> truth_hits;

  // loop over all clusters...
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter) {
    unsigned int cluster_id = *iter;  
    SvtxCluster* cluster = _clustermap->get(cluster_id);

    std::set<PHG4Hit*> new_hits = _clustereval.all_truth_hits(cluster);

    std::set<PHG4Hit*> union_hits; // placeholder for union of new hits and truth hits

    std::set_union(truth_hits.begin(),truth_hits.end(),
		   new_hits.begin(),new_hits.end(),
		   std::inserter(union_hits,union_hits.begin()));

    std::swap(truth_hits,union_hits); // swap union into truth_hits
  }

  if (_do_cache) _cache_all_truth_hits.insert(make_pair(track,truth_hits));
  
  return truth_hits;
}
  
std::set<PHG4Particle*> SvtxTrackEval::all_truth_particles(SvtxTrack* track) {

  if (_do_cache) {
    std::map<SvtxTrack*,std::set<PHG4Particle*> >::iterator iter =
      _cache_all_truth_particles.find(track);
    if (iter != _cache_all_truth_particles.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Particle*> truth_particles;

  // loop over all clusters...
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter) {
    unsigned int cluster_id = *iter;
    SvtxCluster* cluster = _clustermap->get(cluster_id);

    std::set<PHG4Particle*> new_particles = _clustereval.all_truth_particles(cluster);

    std::set<PHG4Particle*> union_particles; // placeholder for union of new particles and truth particles

    std::set_union(truth_particles.begin(),truth_particles.end(),
		   new_particles.begin(),new_particles.end(),
		   std::inserter(union_particles,union_particles.begin()));

    std::swap(truth_particles,union_particles); // swap union into truth_particles
  }

  if (_do_cache) _cache_all_truth_particles.insert(make_pair(track,truth_particles));

  return truth_particles;
}

PHG4Particle* SvtxTrackEval::max_truth_particle_by_nclusters(SvtxTrack* track) {

  if (_do_cache) {
    std::map<SvtxTrack*,PHG4Particle*>::iterator iter =
      _cache_max_truth_particle_by_nclusters.find(track);
    if (iter != _cache_max_truth_particle_by_nclusters.end()) {
      return iter->second;
    }
  }
  
  std::set<PHG4Particle*> particles = all_truth_particles(track);

  PHG4Particle* max_particle = NULL;
  unsigned int max_nclusters = 0;
  
  for(std::set<PHG4Particle*>::iterator iter = particles.begin();
      iter != particles.end();
      ++iter) {
    PHG4Particle* candidate = *iter;
    unsigned int nclusters = get_nclusters_contribution(track,candidate);
    if (nclusters > max_nclusters) {
      max_nclusters = nclusters;
      max_particle = candidate;
    }
  }

  if (_do_cache) _cache_max_truth_particle_by_nclusters.insert(make_pair(track,max_particle));
  
  return max_particle;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Particle* truthparticle) { 

  if (_do_cache) {
    std::map<PHG4Particle*,std::set<SvtxTrack*> >::iterator iter =
      _cache_all_tracks_from_particle.find(truthparticle);
    if (iter !=	_cache_all_tracks_from_particle.end()) {
      return iter->second;
    }
  }
  
  std::set<SvtxTrack*> tracks;
  
  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter) {
    SvtxTrack* track = iter->second;
    
    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
	 iter != track->end_clusters();
	 ++iter) {
      unsigned int cluster_id = *iter;
      SvtxCluster* cluster = _clustermap->get(cluster_id);

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

  if (_do_cache) _cache_all_tracks_from_particle.insert(make_pair(truthparticle,tracks));
  
  return tracks;
}

std::set<SvtxTrack*> SvtxTrackEval::all_tracks_from(PHG4Hit* truthhit) {

  if (_do_cache) {
    std::map<PHG4Hit*,std::set<SvtxTrack*> >::iterator iter =
      _cache_all_tracks_from_g4hit.find(truthhit);
    if (iter != _cache_all_tracks_from_g4hit.end()) {
      return iter->second;
    }
  }
  
  std::set<SvtxTrack*> tracks;
  
  // loop over all SvtxTracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin();
       iter != _trackmap->end();
       ++iter) {
    SvtxTrack* track = iter->second;
    
    // loop over all clusters
    for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
	 iter != track->end_clusters();
	 ++iter) {
      unsigned int cluster_id = *iter;
      SvtxCluster* cluster = _clustermap->get(cluster_id);

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

  if (_do_cache) _cache_all_tracks_from_g4hit.insert(make_pair(truthhit,tracks));

  return tracks;
}

SvtxTrack* SvtxTrackEval::best_track_from(PHG4Particle* truthparticle) { 

  if (_do_cache) {
    std::map<PHG4Particle*,SvtxTrack*>::iterator iter =
      _cache_best_track_from_particle.find(truthparticle);
    if (iter != _cache_best_track_from_particle.end()) {
      return iter->second;
    }
  }

  SvtxTrack* best_track = NULL;
  unsigned int best_count = 0;
  std::set<SvtxTrack*> tracks = all_tracks_from(truthparticle);
  for (std::set<SvtxTrack*>::iterator iter = tracks.begin();
       iter != tracks.end();
       ++iter) {
    SvtxTrack* track = *iter;
    unsigned int count = get_nclusters_contribution(track,truthparticle);
    if (count > best_count) {
      best_track = track;
      best_count = count;
    }
  }
  
  if (_do_cache) _cache_best_track_from_particle.insert(make_pair(truthparticle,best_track));
  
  return best_track;
}

// overlap calculations
unsigned int SvtxTrackEval::get_nclusters_contribution(SvtxTrack* track, PHG4Particle* particle) {

  if (_do_cache) {
    std::map<std::pair<SvtxTrack*,PHG4Particle*>, unsigned int>::iterator iter =
      _cache_get_nclusters_contribution.find(make_pair(track,particle));
    if (iter !=	_cache_get_nclusters_contribution.end()) {
      return iter->second;
    }
  }
  
  unsigned int nclusters = 0; 

  // loop over all clusters
  for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
       iter != track->end_clusters();
       ++iter) {
    unsigned int cluster_id = *iter;
    SvtxCluster* cluster = _clustermap->get(cluster_id);

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
  
  if (_do_cache) _cache_get_nclusters_contribution.insert(make_pair(make_pair(track,particle),nclusters));
  
  return nclusters;
}

void SvtxTrackEval::get_node_pointers(PHCompositeNode *topNode) {

  // need things off of the DST...
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if (!_trackmap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap" << endl;
    exit(-1);
  }
  
  _clustermap = findNode::getClass<SvtxClusterMap>(topNode,"SvtxClusterMap");
  if (!_clustermap) {
    cerr << PHWHERE << " ERROR: Can't find SvtxClusterMap" << endl;
    exit(-1);
  }

  return;
}
