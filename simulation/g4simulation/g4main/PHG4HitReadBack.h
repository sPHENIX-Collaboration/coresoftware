// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4HITREADBACK_H
#define G4MAIN_PHG4HITREADBACK_H

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;

class PHG4HitReadBack : public SubsysReco
{
 public:
  PHG4HitReadBack(const std::string &name="PHG4HITREADBACK");
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
 protected:
};


#endif
