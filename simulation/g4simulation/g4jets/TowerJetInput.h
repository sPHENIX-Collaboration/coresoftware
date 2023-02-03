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
  void set_towerinfo(bool use_towerinfo)
  {
    m_use_towerinfo = use_towerinfo;
  }

  Jet::SRC get_src() override { return _input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

 private:
  Jet::SRC _input;
  bool m_use_towerinfo = true;
  RawTowerDefs::CalorimeterId geocaloid;
};

#endif
