#include "CaloTowerCalib.h"

#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <stdexcept>  // for runtime_error

//____________________________________________________________________________..
CaloTowerCalib::CaloTowerCalib(const std::string &name)
  : SubsysReco(name)
  , m_dettype(CaloTowerCalib::HCALOUT)
  , m_detector("HCALOUT")
  , m_DETECTOR(TowerInfoContainer::HCAL)
  , m_fieldname("")
  , m_runNumber(-1)
{
  std::cout << "CaloTowerCalib::CaloTowerCalib(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloTowerCalib::~CaloTowerCalib()
{
  std::cout << "CaloTowerCalib::~CaloTowerCalib() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloTowerCalib::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);

  EventHeader *evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");

  if (evtHeader)
  {
    m_runNumber = evtHeader->get_RunNumber();
  }
  else
  {
    m_runNumber = -1;
  }

  if (m_dettype == CaloTowerCalib::CEMC)
  {
    m_detector = "CEMC";
    m_DETECTOR = TowerInfoContainer::EMCAL;

    cdb = CDBInterface::instance();
    if (!m_overrideCalibName)
    {
      m_calibName = "cemc_pi0_twrSlope_v1";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "Femc_datadriven_qm1_correction";
    }
    std::string calibdir = cdb->getUrl(m_calibName);
    if (calibdir[0] == '/')
    {
      cdbttree = new CDBTTree(calibdir.c_str());
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file found" << std::endl;
      exit(1);
    }
  }
  else if (m_dettype == CaloTowerCalib::HCALIN)
  {
    m_detector = "HCALIN";
    m_DETECTOR = TowerInfoContainer::HCAL;

    if (!m_overrideCalibName)
    {
      m_calibName = "ihcal_abscalib_cosmic";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "ihcal_abscalib_mip";
    }
    cdb = CDBInterface::instance();
    std::string calibdir = cdb->getUrl(m_calibName);
    if (calibdir[0] == '/')
    {
      cdbttree = new CDBTTree(calibdir.c_str());
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file found" << std::endl;
      exit(1);
    }
  }
  else if (m_dettype == CaloTowerCalib::HCALOUT)
  {
    m_detector = "HCALOUT";
    m_DETECTOR = TowerInfoContainer::HCAL;

    if (!m_overrideCalibName)
    {
      m_calibName = "ohcal_abscalib_cosmic";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "ohcal_abscalib_mip";
    }
    cdb = CDBInterface::instance();
    std::string calibdir = cdb->getUrl(m_calibName);
    if (calibdir[0] == '/')
    {
      cdbttree = new CDBTTree(calibdir.c_str());
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file found" << std::endl;
      exit(1);
    }
  }
  else if (m_dettype == CaloTowerCalib::ZDC)
  {
    m_detector = "ZDC";
    m_DETECTOR = TowerInfoContainer::ZDC;

    if (!m_overrideCalibName)
    {
      m_calibName = "data_driven_zdc_calib";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "zdc_calib";
    }
    cdb = CDBInterface::instance();
    std::string calibdir = cdb->getUrl(m_calibName);
    if (calibdir[0] == '/')
    {
      cdbttree = new CDBTTree(calibdir.c_str());
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file found" << std::endl;
      exit(1);
    }
  }
  
  else if (m_dettype == CaloTowerCalib::EPD)
  {
    m_detector = "EPD";
    m_DETECTOR = TowerInfoContainer::SEPD;
    if (!m_overrideCalibName)
    {
      m_calibName = "noCalibYet";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "noCalibYet";
    }
    cdb = CDBInterface::instance();
    std::string calibdir = cdb->getUrl(m_calibName);
    if (calibdir[0] == '/')
    {
      cdbttree = new CDBTTree(calibdir.c_str());
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file found" << std::endl;
      exit(1);
    }
  }

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode",
                                                           "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  try
  {
    CreateNodeTree(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloTowerCalib::process_event(PHCompositeNode * /*topNode*/)
{
  unsigned int ntowers = _raw_towers->size();
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    unsigned int key = _raw_towers->encode_key(channel);
    TowerInfo *caloinfo_raw = _raw_towers->get_tower_at_channel(channel);
    float raw_amplitude = caloinfo_raw->get_energy();
    short raw_time = caloinfo_raw->get_time();
    float calibconst = cdbttree->GetFloatValue(key, m_fieldname);
    _calib_towers->get_tower_at_channel(channel)->set_energy(raw_amplitude * calibconst);
    if (calibconst == 0){
      _calib_towers->get_tower_at_channel(channel)->set_time(-10);
    }
    else {
      _calib_towers->get_tower_at_channel(channel)->set_time(raw_time);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerCalib::CreateNodeTree(PHCompositeNode *topNode)
{
  std::cout << "creating node" << std::endl;
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerCalibration::CreateNodes");
  }

  // towers
  std::string RawTowerNodeName = "TOWERS_" + m_detector;
  _raw_towers = findNode::getClass<TowerInfoContainerv1>(dstNode,
                                                         RawTowerNodeName);
  if (!_raw_towers)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + RawTowerNodeName + " node in RawTowerCalibration::CreateNodes");
  }

  std::string CalibTowerNodeName = "TOWERINFO_CALIB_" + m_detector;
  _calib_towers = findNode::getClass<TowerInfoContainerv1>(dstNode,
                                                           CalibTowerNodeName);
  if (!_calib_towers)
  {
    _calib_towers = new TowerInfoContainerv1(m_DETECTOR);

    PHIODataNode<PHObject> *calibtowerNode = new PHIODataNode<PHObject>(_calib_towers, CalibTowerNodeName, "PHObject");
    std::cout << "adding calib tower" << std::endl;
    dstNode->addNode(calibtowerNode);
  }

  return;
}
