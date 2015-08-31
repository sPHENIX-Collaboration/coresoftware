
#ifndef __JETTRUTHEVAL_H__
#define __JETTRUTHEVAL_H__

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>

#include <set>
#include <map>

class JetTruthEval {

public:

  JetTruthEval(PHCompositeNode *topNode);
  virtual ~JetTruthEval() {}

  void next_event(PHCompositeNode *topNode);
  void do_caching(bool do_cache) {_do_cache = do_cache;}

  std::set<PHG4Particle*> all_truth_particles(Jet* truthjet);
  std::set<PHG4Hit*>      all_truth_hits(Jet* truthjet);
  
private:

  void get_node_pointers(PHCompositeNode* topNode);
  
  PHG4TruthInfoContainer* _truthinfo;

  bool                                        _do_cache;
  std::set<PHG4Hit*>                          _cache_all_truth_hits;
  std::map<PHG4Particle*,std::set<PHG4Hit*> > _cache_all_truth_hits_g4particle;
  std::map<PHG4Particle*,PHG4Hit*>            _cache_get_innermost_truth_hit;
  std::map<PHG4Particle*,PHG4Hit*>            _cache_get_outermost_truth_hit;
};

#endif // __JETTRUTHEVAL_H__
