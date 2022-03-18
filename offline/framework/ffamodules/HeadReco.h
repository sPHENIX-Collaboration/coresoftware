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
  ~HeadReco() override {}
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 protected:
};

#endif
