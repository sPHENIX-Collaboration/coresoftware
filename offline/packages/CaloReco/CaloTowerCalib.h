// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERCALIB_H
#define CALOTOWERCALIB_H

#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <cassert>
#include <cmath>  // for NAN
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>
#include <string>

class PHCompositeNode;

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
    EPD = 3
  };

  void set_detector_type(CaloTowerCalib::DetectorSystem dettype)
  {
    m_dettype = dettype;
    return;
  }
  void setCalibName(std::string name) 
  {
    m_calibName = name;
    m_overrideCalibName = 1;
    return;
   } 
  void setFieldName(std::string name) 
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
