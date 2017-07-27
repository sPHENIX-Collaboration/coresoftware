#ifndef SVTXTRACKEVAL_H__
#define SVTXTRACKEVAL_H__

#include "SvtxClusterEval.h"
/*
#include "SvtxHitEval.h"
#include "SvtxTruthEval.h"

#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
*/

#include <set>
#include <map>

class PHCompositeNode;

class PHG4Hit;
class PHG4Particle;

class SvtxCluster;
class SvtxClusterMap;
class SvtxTrack;
class SvtxTrackMap;

class SvtxTrackEval {

public:

  SvtxTrackEval(PHCompositeNode *topNode);
  virtual ~SvtxTrackEval();

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {
    _do_cache = do_cache;
    _clustereval.do_caching(do_cache);
  }
  void set_strict(bool strict) {
    _strict = strict;
    _clustereval.set_strict(strict);
  }
  void set_verbosity(int verbosity) {
    _verbosity = verbosity;
    _clustereval.set_verbosity(verbosity);
  }
  
  // access the clustereval (and its cached values)
  SvtxClusterEval* get_cluster_eval() {return &_clustereval;}
  SvtxHitEval*     get_hit_eval() {return _clustereval.get_hit_eval();}
  SvtxTruthEval* get_truth_eval() {return _clustereval.get_truth_eval();}
  
  // backtrace through to PHG4Hits
  std::set<PHG4Hit*>  all_truth_hits (SvtxTrack* track);
  
  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles             (SvtxTrack* track);
  PHG4Particle*           max_truth_particle_by_nclusters (SvtxTrack* track);

  // forwardtrace through to SvtxTracks
  std::set<SvtxTrack*> all_tracks_from(PHG4Particle* truthparticle);
  SvtxTrack*           best_track_from(PHG4Particle* truthparticle);
  std::set<SvtxTrack*> all_tracks_from(PHG4Hit* truthhit);
  std::set<SvtxTrack*> all_tracks_from(SvtxCluster* cluster);
  SvtxTrack*           best_track_from(SvtxCluster* cluster);
  void                 create_cache_track_from_cluster();
  
  // overlap calculations
  void         calc_cluster_contribution(SvtxTrack* svtxtrack, PHG4Particle* truthparticle);  
  unsigned int get_nclusters_contribution(SvtxTrack* svtxtrack, PHG4Particle* truthparticle);  
  unsigned int get_layer_range_contribution(SvtxTrack* track, PHG4Particle* particle, unsigned int start_layer, unsigned int end_layer);
  unsigned int get_nclusters_contribution_by_layer(SvtxTrack* svtxtrack, PHG4Particle* truthparticle);  
  unsigned int get_nwrongclusters_contribution(SvtxTrack* svtxtrack, PHG4Particle* truthparticle);  
  unsigned int get_errors() {return _errors + _clustereval.get_errors();}

private:

  void get_node_pointers(PHCompositeNode* topNode);
  bool has_node_pointers();
  
  SvtxClusterEval _clustereval;
  SvtxTrackMap* _trackmap;
  SvtxClusterMap* _clustermap;

  bool _strict;
  int _verbosity;
  unsigned int _errors;
  
  bool                                                        _do_cache;
  bool                                                        _cache_track_from_cluster_exists;
  std::map<SvtxTrack*,std::set<PHG4Hit*> >                    _cache_all_truth_hits;
  std::map<SvtxTrack*,std::set<PHG4Particle*> >               _cache_all_truth_particles;
  std::map<SvtxTrack*,PHG4Particle*>                          _cache_max_truth_particle_by_nclusters;
  std::map<PHG4Particle*,std::set<SvtxTrack*> >               _cache_all_tracks_from_particle;
  std::map<PHG4Particle*,SvtxTrack* >                         _cache_best_track_from_particle;
  std::map<PHG4Hit*,std::set<SvtxTrack*> >                    _cache_all_tracks_from_g4hit;
  std::map<SvtxCluster*,std::set<SvtxTrack*> >                _cache_all_tracks_from_cluster;
  std::map<SvtxCluster*,SvtxTrack*>                           _cache_best_track_from_cluster;
  std::map<std::pair<SvtxTrack*,PHG4Particle*>, unsigned int> _cache_get_nclusters_contribution;
  std::map<std::pair<SvtxTrack*,PHG4Particle*>, unsigned int> _cache_get_nclusters_contribution_by_layer;
  std::map<std::pair<SvtxTrack*,PHG4Particle*>, unsigned int> _cache_get_nwrongclusters_contribution;
};

#endif // SVTXTRACKEVAL_H__
