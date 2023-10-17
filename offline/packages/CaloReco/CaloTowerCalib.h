// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERCALIB_H
#define CALOTOWERCALIB_H

#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer, TowerIn...

#include <fun4all/SubsysReco.h>

#include <cassert>
#include <iostream>
#include <string>

class CDBInterface;
class CDBTTree;
class PHCompositeNode;
class TowerInfoContainerv1;

class CaloTowerCalib : public SubsysReco
{
 public:
  CaloTowerCalib(const std::string &name = "CaloTowerCalib");

  ~CaloTowerCalib() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void CreateNodeTree(PHCompositeNode *topNode);

  enum DetectorSystem
  {
    CEMC = 0,
    HCALIN = 1,
    HCALOUT = 2,
    EPD = 3,
    ZDC = 4
  };

  void set_detector_type(CaloTowerCalib::DetectorSystem dettype)
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

 private:
  TowerInfoContainerv1 *_raw_towers = nullptr;
  TowerInfoContainerv1 *_calib_towers = nullptr;

  CaloTowerCalib::DetectorSystem m_dettype;

  std::string m_detector;
  TowerInfoContainer::DETECTOR m_DETECTOR;
  std::string m_fieldname;
  std::string m_calibName;
  bool m_overrideCalibName = 0;
  bool m_overrideFieldName = 0;

  CDBInterface *cdb = nullptr;
  CDBTTree *cdbttree = nullptr;
  int m_runNumber;
};

#endif  // CALOTOWERBUILDER_H
