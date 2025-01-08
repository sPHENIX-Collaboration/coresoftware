// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4CALO_G4HITSHIFT_H
#define G4CALO_G4HITSHIFT_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class g4hitshift : public SubsysReco
{
 public:
  g4hitshift(const std::string &name = "g4hitshift");

  ~g4hitshift() override{};

  int process_event(PHCompositeNode *topNode) override;

 private:
};

#endif  // G4CALO_G4HITSHIFT_H
