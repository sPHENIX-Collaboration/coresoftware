
#ifndef __SVTXCLUSTEREVAL_H__
#define __SVTXCLUSTEREVAL_H__

#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxCluster.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>

#include <set>
#include <map>

class SvtxClusterEval {

public:

  SvtxClusterEval(PHCompositeNode *topNode);
  virtual ~SvtxClusterEval() {}

  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits          (SvtxCluster* cluster);
  PHG4Hit*           max_truth_hit_by_energy (SvtxCluster* cluster);
  
  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles          (SvtxCluster* cluster);
  PHG4Particle*           max_truth_particle_by_energy (SvtxCluster* cluster);

  // forwardtrace through to SvtxClusters
  std::set<SvtxCluster*> all_clusters_from(PHG4Particle* truthparticle);
  std::set<SvtxCluster*> all_clusters_from(PHG4Hit* truthhit);
  
  // overlap calculations
  float get_energy_contribution (SvtxCluster* svtxcluster, PHG4Particle* truthparticle);
  
private:
  PHCompositeNode* _topNode;

  std::map<SvtxCluster*,std::set<PHG4Hit*> >            _cache_all_truth_hits;
  std::map<SvtxCluster*,PHG4Hit*>                       _cache_max_truth_hit_by_energy;
  std::map<SvtxCluster*,std::set<PHG4Particle*> >       _cache_all_truth_particles;
  std::map<SvtxCluster*,PHG4Particle* >                 _cache_max_truth_particle_by_energy;
  std::map<PHG4Particle*,std::set<SvtxCluster*> >       _cache_all_clusters_from_particle;
  std::map<PHG4Hit*,std::set<SvtxCluster*> >            _cache_all_clusters_from_g4hit;
  std::map<std::pair<SvtxCluster*,PHG4Particle*>,float> _cache_get_energy_contribution;
};

#endif // __SVTXCLUSTEREVAL_H__
