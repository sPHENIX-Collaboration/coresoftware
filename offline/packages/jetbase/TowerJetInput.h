#ifndef G4JET_TOWERJETINPUT_H
#define G4JET_TOWERJETINPUT_H

#include "JetInput.h"

#include "Jet.h"

#include <iostream>  // for cout, ostream
#include <vector>
#include <calobase/RawTowerDefs.h>
// forward declarations
class PHCompositeNode;
class TowerJetInput : public JetInput
{
 public:
  TowerJetInput(Jet::SRC input);
  ~TowerJetInput() override {}

  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return _input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

 private:
  Jet::SRC _input;
  RawTowerDefs::CalorimeterId geocaloid;
  bool m_use_towerinfo = false;
};

#endif
