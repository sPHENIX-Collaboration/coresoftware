// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAMODULES_TIMING_H
#define FFAMODULES_TIMING_H

#include <fun4all/SubsysReco.h>

#include <string>  // for string

class PHCompositeNode;

class Timing : public SubsysReco
{
 public:
  Timing(const std::string &name = "Timing");
  ~Timing() override {}
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void SetCallCounter(unsigned int i) { calls = i;}
 private:

  unsigned int call_counter {0};
  unsigned int calls {10000};
  unsigned int counter {0};
  time_t starttime {0};
};

#endif
