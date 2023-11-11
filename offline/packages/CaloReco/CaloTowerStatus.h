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
  void setCalibName(const std::string &name)
  {
    m_calibName = name;
    m_overrideCalibName = 1;
    return;
  }
  void setFieldName(const std::string &name)
  {
    m_fieldname = name;
    m_overrideFieldName = 1;
    return;
  }
  void set_inputNodePrefix(const std::string &name)
  {
    m_inputNodePrefix = name;
    return;
  }


 private:
  CaloTowerDefs::DetectorSystem m_dettype;

  std::string m_detector;
  TowerInfoContainer::DETECTOR m_DETECTOR;
  std::string m_fieldname;
  std::string m_calibName;
  bool m_overrideCalibName {false};
  bool m_overrideFieldName {false};
  std::string m_inputNodePrefix {"TOWERS_"};
  std::string RawTowerNodeName;

  TowerInfoContainer *_raw_towers;

  bool m_doHot = 1;

  CDBTTree *cdbttree = nullptr;
  int m_runNumber;

};

#endif  // CALOTOWERBUILDER_H
