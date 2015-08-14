
#ifndef __SVTXTRACKEVAL_H__
#define __SVTXTRACKEVAL_H__

#include "SvtxClusterEval.h"
#include "SvtxHitEval.h"

#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxTrack.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>

#include <set>
#include <map>

class SvtxTrackEval {

public:

  SvtxTrackEval(PHCompositeNode *topNode);
  virtual ~SvtxTrackEval() {}

  void next_event(PHCompositeNode *topNode);
  
  // access the clustereval (and its cached values)
  SvtxClusterEval* get_cluster_eval() {return &_clustereval;}
  SvtxHitEval*     get_hit_eval() {return _clustereval.get_hit_eval();}
  
  // backtrace through to PHG4Hits
  std::set<PHG4Hit*>  all_truth_hits (SvtxTrack* track);
  
  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles             (SvtxTrack* track);
  PHG4Particle*           max_truth_particle_by_nclusters (SvtxTrack* track);

  // forwardtrace through to SvtxTracks
  std::set<SvtxTrack*> all_tracks_from(PHG4Particle* truthparticle);
  std::set<SvtxTrack*> all_tracks_from(PHG4Hit* truthhit);
  
  // overlap calculations
  unsigned int get_nclusters_contribution(SvtxTrack* svtxtrack, PHG4Particle* truthparticle);  
  
private:
  PHCompositeNode* _topNode;
  SvtxClusterEval _clustereval;

  std::map<SvtxTrack*,std::set<PHG4Hit*> >                    _cache_all_truth_hits;
  std::map<SvtxTrack*,std::set<PHG4Particle*> >               _cache_all_truth_particles;
  std::map<SvtxTrack*,PHG4Particle*>                          _cache_max_truth_particle_by_nclusters;
  std::map<PHG4Particle*,std::set<SvtxTrack*> >               _cache_all_tracks_from_particle;
  std::map<PHG4Hit*,std::set<SvtxTrack*> >                    _cache_all_tracks_from_g4hit;
  std::map<std::pair<SvtxTrack*,PHG4Particle*>, unsigned int> _cache_get_nclusters_contribution;
};

#endif // __SVTXTRACKEVAL_H__
