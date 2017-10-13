#ifndef HEPMCFLOWAFTERBURNER_H__
#define HEPMCFLOWAFTERBURNER_H__

#include <fun4all/SubsysReco.h>

#include<string>


//! this class does not seem to be used anymore. One need some work to revive it.
class HepMCFlowAfterBurner : public SubsysReco
{
 public:
 HepMCFlowAfterBurner(const std::string &name = "HEPMCFLOW"): SubsysReco(name) {}
  virtual ~HepMCFlowAfterBurner() {}


 protected:


  int seedset;
  long seed;
  long randomSeed;

};

#endif
