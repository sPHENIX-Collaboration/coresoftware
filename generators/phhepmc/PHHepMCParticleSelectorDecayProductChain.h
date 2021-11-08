#ifndef PHHEPMC_PHHEPMCPARTICLESELECTORDECAYPRODUCTCHAIN_H
#define PHHEPMC_PHHEPMCPARTICLESELECTORDECAYPRODUCTCHAIN_H

#include <fun4all/SubsysReco.h>

#include <string>                // for string
#include <vector>

class PHCompositeNode;

namespace HepMC
{
class GenParticle;
class GenEvent;
}  // namespace HepMC
/// Particle selector for HepMC based events
/// Will write out only _theParticle and _theDaughters (if specified)
/// Special case:  when _theParticle=0, all particles in _theDaughers list
/// will be written out no matter where they come from
class PHHepMCParticleSelectorDecayProductChain : public SubsysReco
{
 public:
  PHHepMCParticleSelectorDecayProductChain(const std::string& name = "PARTICLESELECTOR");
  ~PHHepMCParticleSelectorDecayProductChain() override {}

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  /// Set the ID of the particle you want in your output.
  virtual void SetParticle(const int pid);

  /// Add an ancestor of the particle you want in your output.
  virtual void AddAncestor(const int pid);

  /// Add decay products of the particle you want in your output.
  virtual void AddDaughter(const int pid);

  //! embedding ID for the event to be processed
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return _embedding_id; }
  //
  //! embedding ID for the event to be processed
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { _embedding_id = id; }

 protected:
  /// find out if a particle comes from one of _theAncestors
  HepMC::GenParticle* GetParent(HepMC::GenParticle* p, HepMC::GenEvent* event);

  /// The particle you want to have in your output
  int _theParticle;
  /// List of possible decay products of the particle you want in your output
  /// Ignored if empty
  std::vector<int> _theDaughters;
  /// List of possible ancestors of the particle you want in your output
  /// Ignored if empty
  std::vector<int> _theAncestors;

  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id;
};

#endif
