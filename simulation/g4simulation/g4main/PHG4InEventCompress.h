// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4INEVENTCOMPRESS_H
#define G4MAIN_PHG4INEVENTCOMPRESS_H

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;
class VariableArray;

class PHG4InEventCompress: public SubsysReco
{
 public:
  PHG4InEventCompress(const std::string &name = "PHG4InEventCompress");
  ~PHG4InEventCompress() override {}
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 protected:
  VariableArray *vtxarray;
  VariableArray *particlearray;
};

#endif
