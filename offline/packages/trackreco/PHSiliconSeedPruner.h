// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef PHSILICONSEEDPRUNER_H
#define PHSILICONSEEDPRUNER_H

#include <fun4all/SubsysReco.h>

#include <string>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>

class PHCompositeNode;

class PHSiliconSeedPruner : public SubsysReco
{
 public:
  PHSiliconSeedPruner(const std::string& name = "PHSiliconSeedPruner");

  ~PHSiliconSeedPruner() override;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;
};

#endif  // PHSILICONSEEDPRUNER_H
