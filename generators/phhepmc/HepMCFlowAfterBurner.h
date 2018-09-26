#ifndef PHHEPMC_HEPMCFLOWAFTERBURNER_H
#define PHHEPMC_HEPMCFLOWAFTERBURNER_H

#include <fun4all/SubsysReco.h>

#include<string>

class HepMCFlowAfterBurner : public SubsysReco
{
 public:
  HepMCFlowAfterBurner(const std::string &name = "HEPMCFLOW");
  virtual ~HepMCFlowAfterBurner() {}

  int Init(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  void setConfigFileName(const std::string &fnam)
  {
    config_filename = fnam;
  }
  void setAlgorithmName(const std::string &name)
  {
    algorithmName_macro = name;
    algorithmName = algorithmName_macro;
  }
  void setMinEta(const float f)
  {
    mineta_macro = f;
    mineta = mineta_macro;
  }
  void setMaxEta(const float f)
  {
    maxeta_macro = f;
    maxeta = maxeta_macro;
  }
  void setMinPt(const float f)
  {
    minpt_macro = f;
    minpt = minpt_macro;
  }
  void setMaxPt(const float f)
  {
    maxpt_macro = f;
    maxpt = maxpt_macro;
  }
  void setSeed(const long il);

  void SaveRandomState(const std::string &savefile = "HepMCFlowAfterBurner.ransave");
  void RestoreRandomState(const std::string &savefile = "HepMCFlowAfterBurner.ransave");

 protected:

  std::string config_filename;
  std::string algorithmName;
  std::string algorithmName_macro;

  float mineta;
  float maxeta;

  float mineta_macro;
  float maxeta_macro;

  float minpt;
  float maxpt;

  float minpt_macro;
  float maxpt_macro;

  int seedset;
  long seed;
  long randomSeed;


};

#endif
