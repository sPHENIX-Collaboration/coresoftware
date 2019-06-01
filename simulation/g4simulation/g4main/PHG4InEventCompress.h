// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4INEVENTCOMPRESS_H
#define G4MAIN_PHG4INEVENTCOMPRESS_H

#include <fun4all/SubsysReco.h>

class VariableArray;

class PHG4InEventCompress: public SubsysReco
{
 public:
  PHG4InEventCompress(const std::string &name = "PHG4InEventCompress");
  virtual ~PHG4InEventCompress() {}
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 protected:
  VariableArray *vtxarray;
  VariableArray *particlearray;
};

#endif
