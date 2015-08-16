
#ifndef __SVTXHITEVAL_H__
#define __SVTXHITEVAL_H__

#include <phool/PHCompositeNode.h>
#include <g4hough/SvtxHit.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>

#include <set>
#include <map>

class SvtxHitEval {

public:

  SvtxHitEval(PHCompositeNode *topNode);
  virtual ~SvtxHitEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}

  PHG4CylinderCell* get_cell(SvtxHit* hit);
  
  // backtrace through to PHG4Hits
  std::set<PHG4Hit*> all_truth_hits          (SvtxHit* hit);
  PHG4Hit*           max_truth_hit_by_energy (SvtxHit* hit);
  
  // backtrace through to PHG4Particles
  std::set<PHG4Particle*> all_truth_particles          (SvtxHit* hit);
  PHG4Particle*           max_truth_particle_by_energy (SvtxHit* hit);

  // forwardtrace through to SvtxHits
  std::set<SvtxHit*> all_hits_from(PHG4Particle* truthparticle);
  std::set<SvtxHit*> all_hits_from(PHG4Hit* truthhit);
  SvtxHit*           best_hit_from(PHG4Hit* truthhit);
  
  // overlap calculations
  float get_energy_contribution (SvtxHit* svtxhit, PHG4Particle* truthparticle);
  float get_energy_contribution (SvtxHit* svtxhit, PHG4Hit* truthhit);
  
private:
  PHCompositeNode* _topNode;

  bool                                              _do_cache;
  std::map<SvtxHit*,std::set<PHG4Hit*> >            _cache_all_truth_hits;
  std::map<SvtxHit*,PHG4Hit*>                       _cache_max_truth_hit_by_energy;
  std::map<SvtxHit*,std::set<PHG4Particle*> >       _cache_all_truth_particles;
  std::map<SvtxHit*,PHG4Particle* >                 _cache_max_truth_particle_by_energy;
  std::map<PHG4Particle*,std::set<SvtxHit*> >       _cache_all_hits_from_particle;
  std::map<PHG4Hit*,std::set<SvtxHit*> >            _cache_all_hits_from_g4hit;
  std::map<PHG4Hit*,SvtxHit*>                       _cache_best_hit_from_g4hit;
  std::map<std::pair<SvtxHit*,PHG4Particle*>,float> _cache_get_energy_contribution_g4particle;
  std::map<std::pair<SvtxHit*,PHG4Hit*>,float>      _cache_get_energy_contribution_g4hit;
};

#endif // __SVTXHITEVAL_H__
