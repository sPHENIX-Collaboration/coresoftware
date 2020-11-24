// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAMODULES_HEADRECO_H
#define FFAMODULES_HEADRECO_H

#include <fun4all/SubsysReco.h>

#include <string>  // for string

class PHCompositeNode;

class HeadReco : public SubsysReco
{
 public:
  HeadReco(const std::string &name = "HeadReco");
  virtual ~HeadReco() {}
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int EndRun(const int runno);

 protected:
};

#endif
