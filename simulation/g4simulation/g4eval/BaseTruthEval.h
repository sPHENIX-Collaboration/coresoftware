#ifndef G4EVAL_BASETRUTHEVAL_H
#define G4EVAL_BASETRUTHEVAL_H

#include <set>

class PHCompositeNode;
class PHG4Hit;
class PHG4Particle;
class PHG4Shower;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;

class BaseTruthEval
{
 public:
  explicit BaseTruthEval(PHCompositeNode* topNode);
  virtual ~BaseTruthEval();

  /// reinitialize the eval for a new event
  void next_event(PHCompositeNode* topNode);

  /// strict mode will assert when an error is detected
  /// non-strict mode will notice and report at the End()
  void set_strict(bool strict) { m_Strict = strict; }

  /// get a count of the errors discovered thus far
  unsigned int get_errors() { return m_Errors; }

  /// adjust the messaging from the evalutaion module
  void set_verbosity(int verbosity) { m_Verbosity = verbosity; }

  // ---reduced sim node or better---------------------------------------------

  bool has_reduced_node_pointers();

  /// what was the embed flag set for this particle?
  int get_embed(PHG4Particle* particle);

  /// what was the vertex creation point of the particle?
  PHG4VtxPoint* get_vertex(PHG4Particle* particle);

  /// is this a primary shower?
  bool is_primary(PHG4Shower* shower);

  /// is this a primary particle?
  bool is_primary(PHG4Particle* particle);

  /// what was the primary shower for this possibly secondary shower?
  PHG4Shower* get_primary_shower(PHG4Shower* shower);

  /// what was the primary shower that is associated with this particle?
  PHG4Shower* get_primary_shower(PHG4Particle* particle);

  /// what was the primary particle that is associated with this particle?
  PHG4Particle* get_primary_particle(PHG4Particle* particle);

  /// what was the primary particle that is associated with this shower?
  PHG4Particle* get_primary_particle(PHG4Shower* shower);

  /// which secondary showers are inside this shower?
  std::set<PHG4Shower*> all_secondary_showers(PHG4Shower* shower);

  /// do these two shower pointers resolve to the same shower?
  bool are_same_shower(PHG4Shower* s1, PHG4Shower* s2);

  /// do these two particle pointers resolve to the same particle?
  bool are_same_particle(PHG4Particle* p1, PHG4Particle* p2);

  /// do these two vertex pointers resolve to the same vertex?
  bool are_same_vertex(PHG4VtxPoint* vtx1, PHG4VtxPoint* vtx2);

  // ---full sim node required--------------------------------------------------

  bool has_full_node_pointers() { return has_reduced_node_pointers(); }

  /// which particle left this truth hit?
  PHG4Particle* get_particle(PHG4Hit* g4hit);

  /// which primary shower contains this truth hit?
  PHG4Shower* get_primary_shower(PHG4Hit* g4hit);

  /// which primary particle resulted in this truth hit?
  PHG4Particle* get_primary_particle(PHG4Hit* g4hit);

  /// which is the particle associated with this track ID?
  PHG4Particle* get_particle(const int trackid); 

  /// is this truth hit inside this shower?
  bool is_g4hit_from_primary_shower(PHG4Hit* g4hit, PHG4Shower* shower);

  /// was this truth hit left by this particle?
  bool is_g4hit_from_particle(PHG4Hit* g4hit, PHG4Particle* particle);

 private:
  void get_node_pointers(PHCompositeNode* topNode);

  PHG4TruthInfoContainer* m_TruthInfo;

  bool m_Strict;
  int m_Verbosity;
  unsigned int m_Errors;
};

#endif  // G4EVAL_BASETRUTHEVAL_H
