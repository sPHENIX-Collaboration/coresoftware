
#ifndef __CALOTRUTHEVAL_H__
#define __CALOTRUTHEVAL_H__

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <string>
#include <set>
#include <map>

class CaloTruthEval {

public:

  CaloTruthEval(PHCompositeNode *topNode, std::string caloname = "CEMC"); // CEMC, HCALIN, HCALOUT
  virtual ~CaloTruthEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}

  std::set<PHG4Hit*> all_truth_hits(PHG4Particle* particle);  
  PHG4Particle*      get_parent_particle(PHG4Hit* g4hit);
  PHG4Particle*      get_primary_particle(PHG4Hit* g4hit);  
  int                get_embed(PHG4Particle* particle);
  PHG4VtxPoint*      get_vertex(PHG4Particle* particle);

  bool               is_primary(PHG4Particle* particle);
  std::set<PHG4Hit*> get_shower_from_primary(PHG4Particle* primary);  
  float              get_shower_moliere_radius(PHG4Particle* primary);
  float              get_shower_energy_deposit(PHG4Particle* primary);
  
private:
  PHCompositeNode* _topNode;
  std::string _caloname;

  bool                                        _do_cache;
  std::map<PHG4Particle*,std::set<PHG4Hit*> > _cache_all_truth_hits_g4particle;
  std::map<PHG4Particle*,float>               _cache_get_shower_moliere_radius;
  std::map<PHG4Particle*,float>               _cache_get_shower_sample_fraction;
};

#endif // __CALOTRUTHEVAL_H__
