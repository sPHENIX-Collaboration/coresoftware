#ifndef G4EVAL_CALOTRUTHEVAL_H
#define G4EVAL_CALOTRUTHEVAL_H

#include "BaseTruthEval.h"

#include <map>
#include <set>
#include <string>

class PHCompositeNode;

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4Shower;
class PHG4VtxPoint;

class CaloTruthEval
{
 public:
  /// example caloname: CEMC, HCALIN, HCALOUT
  CaloTruthEval(PHCompositeNode* topNode, const std::string& caloname);
  virtual ~CaloTruthEval();

  /// get the hash id for this calorimeter volume
  int get_caloid() { return _caloid; }

  /// reinitialize the eval for a new event
  void next_event(PHCompositeNode* topNode);

  /// activate or deactivate the memory caching inside the evaluation module
  void do_caching(bool do_cache) { _do_cache = do_cache; }

  /// strict mode will assert when an error is detected
  /// non-strict mode will notice and report at the End()
  void set_strict(bool strict)
  {
    _strict = strict;
    _basetrutheval.set_strict(strict);
  }

  /// get a count of the errors discovered thus far
  unsigned int get_errors() { return _errors; }

  /// adjust the messaging from the evalutaion module
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _basetrutheval.set_verbosity(verbosity);
  }

  // ---reduced sim node or better----------------------------------------------

  /// has the eval initialized correctly for reduced sim DST nodes?
  bool has_reduced_node_pointers();

  /// what was the primary shower associated with the shower object?
  PHG4Shower* get_primary_shower(PHG4Shower* shower);

  /// what was the primary shower associated with the shower object?
  PHG4Shower* get_primary_shower(PHG4Particle* particle);

  /// what was the primary particle associated with the shower object?
  PHG4Particle* get_primary_particle(PHG4Shower* shower);

  /// what was the primary particle associated with another or same particle?
  PHG4Particle* get_primary_particle(PHG4Particle* particle);

  /// what was the embed flag passed with this particle?
  int get_embed(PHG4Particle* particle);

  /// what was particle's creation point?
  PHG4VtxPoint* get_vertex(PHG4Particle* particle);

  /// are these two pointers in fact the same shower?
  bool are_same_shower(PHG4Shower* s1, PHG4Shower* s2);

  /// are these two pointers in fact the same particle?
  bool are_same_particle(PHG4Particle* p1, PHG4Particle* p2);

  /// are these two pointers in fact the same vertex?
  bool are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2);

  /// is this a primary shower?
  bool is_primary(PHG4Shower* shower);

  /// is this a primary particle?
  bool is_primary(PHG4Particle* particle);

  /// how much energy did this primary and its shower deposit in the calo volume?
  float get_shower_energy_deposit(PHG4Particle* primary);

  /// what was the electron/hadron radio of energy deposits inside the calo colume?
  float get_shower_eh_ratio(PHG4Particle* primary);

  // ---full sim node required--------------------------------------------------

  /// has the eval initialized correctly for full sim DST nodes?
  bool has_full_node_pointers();

  /// what truth hits are contained within this shower entry in this calo volume?
  std::set<PHG4Hit*> all_truth_hits(PHG4Shower* shower);

  /// what truth hits were left behind by this particle in this calo volume?
  std::set<PHG4Hit*> all_truth_hits(PHG4Particle* particle);

  /// what particle created this truth hit?
  PHG4Particle* get_parent_particle(PHG4Hit* g4hit);

  /// what primary particle was responsible for this truth hit?
  PHG4Particle* get_primary_particle(PHG4Hit* g4hit);

  /// did this particle create this truth hit?
  bool is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle);

  /// what truth hits are the result of this primary particle and its shower
  std::set<PHG4Hit*> get_shower_hits_from_primary(PHG4Particle* primary);

 private:
  void get_node_pointers(PHCompositeNode* topNode);

  BaseTruthEval _basetrutheval;

  std::string _caloname;
  int _caloid;
  PHG4TruthInfoContainer* _truthinfo;
  PHG4HitContainer* _g4hits;
  int _g4hit_container_id;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::map<PHG4Particle*, float> _cache_get_shower_energy_deposit;
  std::map<PHG4Shower*, std::set<PHG4Hit*> > _cache_all_truth_hits_g4shower;
  std::map<PHG4Particle*, std::set<PHG4Hit*> > _cache_all_truth_hits_g4particle;
  std::map<PHG4Hit*, PHG4Particle*> _cache_get_primary_particle_g4hit;
  std::map<PHG4Particle*, std::set<PHG4Hit*> > _cache_get_shower_hits_from_primary;
};

#endif  // G4EVAL_CALOTRUTHEVAL_H
