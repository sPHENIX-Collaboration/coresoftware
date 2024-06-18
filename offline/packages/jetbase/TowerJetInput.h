#ifndef JETBASE_TOWERJETINPUT_H
#define JETBASE_TOWERJETINPUT_H

#include "Jet.h"
#include "JetInput.h"

#include <calobase/RawTowerDefs.h>

#include <iostream>  // for cout, ostream

#include <vector>
// forward declarations
class PHCompositeNode;
class TowerJetInput : public JetInput
{
 public:
  TowerJetInput(Jet::SRC input, const std::string &prefix = "TOWERINFO_CALIB");
  ~TowerJetInput() override {}

  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return m_input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

 private:
  Jet::SRC m_input;
  RawTowerDefs::CalorimeterId geocaloid{RawTowerDefs::CalorimeterId::NONE};
  bool m_use_towerinfo = false;
  std::string m_towerNodePrefix;
  std::string towerName;
};

#endif
