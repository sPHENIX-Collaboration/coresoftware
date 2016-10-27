#ifndef PHHepMCParticleSelectorDecayProductChain_H__
#define PHHepMCParticleSelectorDecayProductChain_H__

#include "PHHepMCParticleSelectorDecayProduct.h"

/// Particle selector for B -> J/spi -> e+e- process
/// Will write out B-mason, J/psi from that B-meson decay and 
/// electrons/positrons from J/psi decay
class PHHepMCParticleSelectorDecayProductChain: public PHHepMCParticleSelectorDecayProduct
{
 public:
  PHHepMCParticleSelectorDecayProductChain(const std::string &name="B2JPSISELECTOR");
  virtual ~PHHepMCParticleSelectorDecayProductChain(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

 protected:

};

#endif


