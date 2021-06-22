// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4INEVENTREADBACK_H
#define G4MAIN_PHG4INEVENTREADBACK_H

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;

class PHG4InEventReadBack: public SubsysReco
{
 public:
  PHG4InEventReadBack(const std::string &name = "PHG4InEventReadBack");
  ~PHG4InEventReadBack() override {}
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 protected:
  
};

#endif
