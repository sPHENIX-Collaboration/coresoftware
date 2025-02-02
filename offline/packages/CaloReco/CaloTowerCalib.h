// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERCALIB_H
#define CALOTOWERCALIB_H

#include "CaloTowerDefs.h"

#include <calobase/TowerInfoContainer.h>  // for TowerInfoContainer, TowerIn...

#include <fun4all/SubsysReco.h>

#include <iostream>
#include <string>

class CDBTTree;
class PHCompositeNode;
class TowerInfoContainer;

class CaloTowerCalib : public SubsysReco
{
 public:
  CaloTowerCalib(const std::string &name = "CaloTowerCalib");

  ~CaloTowerCalib() override;

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
  void set_outputNodePrefix(const std::string &name)
  {
    m_outputNodePrefix = name;
    return;
  }

  void set_directURL(const std::string &url)
  {
    m_giveDirectURL = true;
    m_directURL = url;
  }

  void set_directURL_timecalib(const std::string &url)
  {
    m_giveDirectURL_time = true;
    m_directURL_time = url;
  }

  void set_directURL_ZScrosscalib(const std::string &url)
  {
    m_giveDirectURL_ZScrosscalib = true;
    m_directURL_ZScrosscalib = url;
  }

  void set_doZScrosscalib(bool doZScrosscalib)
  {
    m_doZScrosscalib = doZScrosscalib;
  }

  void set_dotimecalib(bool dotimecalib)
  {
    m_dotimecalib = dotimecalib;
  }

  void set_doCalibOnly(bool docalib = true)
  {
    if (docalib) 
    {
      m_dotimecalib = false;
      m_doZScrosscalib = false;
    }
  }

  void set_use_TowerInfov2(bool use) { m_use_TowerInfov2 = use; }

 private:
  CaloTowerDefs::DetectorSystem m_dettype;

  std::string m_detector;
  TowerInfoContainer::DETECTOR m_DETECTOR;
  std::string m_fieldname;
  std::string m_calibName;
  std::string m_fieldname_time;
  std::string m_calibName_time;
  std::string m_fieldname_ZScrosscalib;
  std::string m_calibName_ZScrosscalib;
  bool m_overrideCalibName{false};
  bool m_overrideFieldName{false};
  std::string m_inputNodePrefix{"TOWERS_"};
  std::string m_outputNodePrefix{"TOWERINFO_CALIB_"};
  std::string RawTowerNodeName;
  std::string CalibTowerNodeName;

  bool m_use_TowerInfov2 = 0;

  bool m_giveDirectURL = false;
  std::string m_directURL = "";

  bool m_giveDirectURL_time = false;
  std::string m_directURL_time = "";
  bool m_dotimecalib = true;

  bool m_giveDirectURL_ZScrosscalib = false;
  std::string m_directURL_ZScrosscalib = "";
  bool m_doZScrosscalib = true;

  CDBTTree *cdbttree = nullptr;
  CDBTTree *cdbttree_time = nullptr;
  CDBTTree *cdbttree_ZScrosscalib = nullptr;
  int m_runNumber;
};

#endif  // CALOTOWERBUILDER_H
