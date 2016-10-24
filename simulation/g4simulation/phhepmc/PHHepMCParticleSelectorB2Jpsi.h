#ifndef PHHepMCParticleSelectorB2Jpsi_H__
#define PHHepMCParticleSelectorB2Jpsi_H__

#include "PHHepMCParticleSelectorBase.h"

class PHHepMCParticleSelectorB2Jpsi: public PHHepMCParticleSelectorBase
{
 public:
  PHHepMCParticleSelectorB2Jpsi(const std::string &name="B2JPSISELECTOR");
  virtual ~PHHepMCParticleSelectorB2Jpsi(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

 protected:

};

#endif


