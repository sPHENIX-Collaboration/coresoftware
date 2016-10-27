#ifndef PHHepMCParticleSelectorB2Jpsi_H__
#define PHHepMCParticleSelectorB2Jpsi_H__

#include "PHHepMCParticleSelectorDecayProduct.h"

/// Particle selector for B -> J/spi -> e+e- process
/// Will write out B-mason, J/psi from that B-meson decay and 
/// electrons/positrons from J/psi decay
class PHHepMCParticleSelectorB2Jpsi: public PHHepMCParticleSelectorDecayProduct
{
 public:
  PHHepMCParticleSelectorB2Jpsi(const std::string &name="B2JPSISELECTOR");
  virtual ~PHHepMCParticleSelectorB2Jpsi(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

 protected:

};

#endif


