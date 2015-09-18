
#ifndef __CALOTRUTHEVAL_H__
#define __CALOTRUTHEVAL_H__

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <string>
#include <set>
#include <map>

class CaloTruthEval {

public:

  CaloTruthEval(PHCompositeNode *topNode, std::string caloname); // CEMC, HCALIN, HCALOUT
  virtual ~CaloTruthEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}

  std::set<PHG4Hit*> all_truth_hits(PHG4Particle* particle);  
  PHG4Particle*      get_parent_particle(PHG4Hit* g4hit);
  PHG4Particle*      get_primary_particle(PHG4Particle* particle);
  PHG4Particle*      get_primary_particle(PHG4Hit* g4hit);  
  int                get_embed(PHG4Particle* particle);
  PHG4VtxPoint*      get_vertex(PHG4Particle* particle);

  bool               is_primary(PHG4Particle* particle);
  std::set<PHG4Hit*> get_shower_from_primary(PHG4Particle* primary);  
  float              get_shower_moliere_radius(PHG4Particle* primary);
  float              get_shower_energy_deposit(PHG4Particle* primary);
  
private:

  void get_node_pointers(PHCompositeNode *topNode);
  
  std::string _caloname;
  PHG4TruthInfoContainer* _truthinfo;
  PHG4HitContainer* _g4hits;

  bool                                        _do_cache;
  std::map<PHG4Particle*,std::set<PHG4Hit*> > _cache_all_truth_hits_g4particle;
  std::map<PHG4Particle*,PHG4Particle*>       _cache_get_primary_particle_g4particle;
  std::map<PHG4Hit*,PHG4Particle*>            _cache_get_primary_particle_g4hit;
  //std::map<PHG4Particle*,bool>                _cache_is_primary;
  std::map<PHG4Particle*,std::set<PHG4Hit*> > _cache_get_shower_from_primary;
  std::map<PHG4Particle*,float>               _cache_get_shower_moliere_radius;
  std::map<PHG4Particle*,float>               _cache_get_shower_energy_deposit;
};

#endif // __CALOTRUTHEVAL_H__
