// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HEADRECO_H
#define G4MAIN_PHG4HEADRECO_H

#include <fun4all/SubsysReco.h>

#include <string>  // for string

class PHCompositeNode;

class PHG4HeadReco : public SubsysReco
{
 public:
  PHG4HeadReco(const std::string &name = "PHG4HeadReco");
  ~PHG4HeadReco() override {}
  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 protected:
  int evtseq = 0;
};

#endif
