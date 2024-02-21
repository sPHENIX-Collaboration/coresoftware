// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERSTATUS_H
#define CALOTOWERSTATUS_H

#include "CaloTowerDefs.h"

#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer, TowerIn...

#include <fun4all/SubsysReco.h>

#include <cassert>
#include <iostream>
#include <string>

class CDBTTree;
class PHCompositeNode;
class TowerInfoContainer;

class CaloTowerStatus : public SubsysReco
{
 public:
  CaloTowerStatus(const std::string &name = "CaloTowerStatus");

  ~CaloTowerStatus() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  void set_detector_type(CaloTowerDefs::DetectorSystem dettype)
  {
    m_dettype = dettype;
    return;
  }
  void set_inputNodePrefix(const std::string &name)
  {
    m_inputNodePrefix = name;
    return;
  }
  void set_badChi2_treshold(float threshold)
  {
    badChi2_treshold = threshold;
    return;
  }
  void set_fraction_badChi2_threshold(float threshold)
  {
    fraction_badChi2_threshold = threshold;
    return;
  }
  void set_time_cut(float threshold)
  {
    time_cut = threshold;
    return;
  }

 private:
  TowerInfoContainer *m_raw_towers{nullptr};

  CDBTTree *m_cdbttree_chi2{nullptr};
  CDBTTree *m_cdbttree_time{nullptr};
  CDBTTree *m_cdbttree_hotMap{nullptr};

  bool m_doHotChi2{true};
  bool m_doTime{true};
  bool m_doHotMap{true};

  CaloTowerDefs::DetectorSystem m_dettype{CaloTowerDefs::DETECTOR_INVALID};

  std::string m_detector;
  std::string m_fieldname_time;
  std::string m_calibName_time;
  std::string m_fieldname_chi2;
  std::string m_calibName_chi2;
  std::string m_fieldname_hotMap;
  std::string m_calibName_hotMap;
  std::string m_inputNodePrefix{"TOWERS_"};

  float badChi2_treshold = 1e4;
  float fraction_badChi2_threshold = 0.01;
  float time_cut = 2; // number of samples from the mean time for the channel in the run
};

#endif  // CALOTOWERBUILDER_H
