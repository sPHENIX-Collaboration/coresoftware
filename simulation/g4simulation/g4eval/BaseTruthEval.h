
#ifndef __BASETRUTHEVAL_H__
#define __BASETRUTHEVAL_H__

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Shower.h>

class BaseTruthEval {

public:

  BaseTruthEval(PHCompositeNode *topNode);
  virtual ~BaseTruthEval();

  void next_event(PHCompositeNode *topNode);
  void set_strict(bool strict) {_strict = strict;}
  void set_verbosity(int verbosity) {_verbosity = verbosity;}
  
  PHG4Particle*      get_particle(PHG4Hit* g4hit);
  int                get_embed(PHG4Particle* particle);
  PHG4VtxPoint*      get_vertex(PHG4Particle* particle);

  bool                  is_primary(PHG4Particle* particle);
  bool                  is_primary(PHG4Shower* shower);
  PHG4Particle*         get_primary(PHG4Hit* g4hit);  
  PHG4Particle*         get_primary(PHG4Particle* particle);  
  PHG4Particle*         get_primary(PHG4Shower* shower);
  PHG4Shower*           get_shower_object_from_primary(PHG4Particle* particle);
  std::set<PHG4Shower*> all_subshower_objects(PHG4Shower* shower);
  
  bool               is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle);
  bool               are_same_particle(PHG4Particle* p1, PHG4Particle* p2);
  bool               are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2);

  unsigned int       get_errors() {return _errors;}
  
private:

  void get_node_pointers(PHCompositeNode* topNode);
  bool has_node_pointers();
  
  PHG4TruthInfoContainer* _truthinfo;

  bool _strict;
  int _verbosity;
  unsigned int _errors;
};

#endif // __BASETRUTHEVAL_H__
