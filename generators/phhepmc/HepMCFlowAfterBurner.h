// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef PHHEPMC_HEPMCFLOWAFTERBURNER_H
#define PHHEPMC_HEPMCFLOWAFTERBURNER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class HepMCFlowAfterBurner : public SubsysReco
{
 public:
  HepMCFlowAfterBurner(const std::string &name = "HEPMCFLOWAFTERBURNER");
  ~HepMCFlowAfterBurner() override {}

  int Init(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;

  void Print(const std::string &what = "ALL") const override;

  void setAlgorithmName(const std::string &name);

  void setMinEta(const float f)
  {
    mineta = f;
  }
  void setMaxEta(const float f)
  {
    maxeta = f;
  }
  void setMinPt(const float f)
  {
    minpt = f;
  }
  void setMaxPt(const float f)
  {
    maxpt = f;
  }
  void setSeed(const long il);

  void SaveRandomState(const std::string &savefile = "HepMCFlowAfterBurner.ransave");
  void RestoreRandomState(const std::string &savefile = "HepMCFlowAfterBurner.ransave");

 protected:
  std::string config_filename;
  std::string algorithmName;

  float mineta;
  float maxeta;

  float minpt;
  float maxpt;

  int seedset;
  long seed;
  long randomSeed;
};

#endif
