// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4CALO_G4HITSHIFTHCAL_H
#define G4CALO_G4HITSHIFTHCAL_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class g4hitshifthcal : public SubsysReco
{
 public:

  g4hitshifthcal(const std::string &name = "g4hitshifthcal");

  ~g4hitshifthcal() override = default;

  int process_event(PHCompositeNode *topNode) override;

 private:
};

#endif // G4CALO_G4HITSHIFTHCAL_H
