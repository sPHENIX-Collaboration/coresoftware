#ifndef PHHepMCParitcleSelectorBase_H_
#define PHHepMCParitcleSelectorBase_H_

#include <fun4all/SubsysReco.h>

#include <vector>

class PHHepMCParticleSelectorBase: public SubsysReco
{
 public:
  virtual ~PHHepMCParticleSelectorBase();

  virtual int InitRun(PHCompositeNode *topNode);
  virtual int process_event(PHCompositeNode *topNode);

  virtual void SetParticle(const int pid);
  virtual void AddParent(const int pid);
  virtual void AddDaughter(const int pid);

  void SetTriggerParticle(int part) { _theTrigger=part; }

 protected:
  PHHepMCParticleSelectorBase(const std::string &name="PARTICLESELECTORBASE");

  int _theParticle;
  std::vector<int> _theDaughters;
  std::vector<int> _theParents;

  int _theTrigger;

};

#endif


