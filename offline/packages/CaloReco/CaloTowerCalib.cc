#include "CaloTowerCalib.h"
#include "CaloTowerDefs.h"

#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfov2.h>

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

#include <TSystem.h>

#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <stdexcept>  // for runtime_error

//____________________________________________________________________________..
CaloTowerCalib::CaloTowerCalib(const std::string &name)
  : SubsysReco(name)
  , m_dettype(CaloTowerDefs::HCALOUT)
  , m_detector("HCALOUT")
  , m_DETECTOR(TowerInfoContainer::HCAL)
  , m_fieldname("")
  , m_runNumber(-1)
{
  if (Verbosity() > 0) std::cout << "CaloTowerCalib::CaloTowerCalib(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloTowerCalib::~CaloTowerCalib()
{
  delete cdbttree;
  if (Verbosity() > 0) std::cout << "CaloTowerCalib::~CaloTowerCalib() Calling dtor" << std::endl;
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

  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_detector = "CEMC";
    m_DETECTOR = TowerInfoContainer::EMCAL;
    std::string default_time_independent_calib = "cemc_pi0_twrSlope_v1_default";

    if (!m_overrideCalibName)
    {
      m_calibName = "cemc_pi0_twrSlope_v1";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "Femc_datadriven_qm1_correction";
    }
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      calibdir = CDBInterface::instance()->getUrl(default_time_independent_calib);

      if (calibdir.empty())
      {
        std::cout << "CaloTowerCalib::::InitRun No EMCal Calibration NOT even a default" << std::endl;
        exit(1);
      }
      cdbttree = new CDBTTree(calibdir);
      std::cout << "CaloTowerCalib::::InitRun No specific file for " << m_calibName << " found, using default calib " << default_time_independent_calib << std::endl;
    }
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
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
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file for domain " << m_calibName << " found" << std::endl;
      exit(1);
    }
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
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
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file for domain " << m_calibName << " found" << std::endl;
      exit(1);
    }
  }
  else if (m_dettype == CaloTowerDefs::ZDC)
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
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file for domain " << m_calibName << " found" << std::endl;
      exit(1);
    }
  }

  else if (m_dettype == CaloTowerDefs::SEPD)
  {
    m_detector = "SEPD";
    m_DETECTOR = TowerInfoContainer::SEPD;
    if (!m_overrideCalibName)
    {
      m_calibName = "noCalibYet";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "noCalibYet";
    }
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerCalib::::InitRun No calibration file for domain " << m_calibName << " found" << std::endl;
      exit(1);
    }
  }

  if (m_giveDirectURL)
  {
    cdbttree = new CDBTTree(m_directURL);
  }

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
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
  if (Verbosity() > 0) topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloTowerCalib::process_event(PHCompositeNode *topNode)
{
  TowerInfoContainer *_raw_towers = findNode::getClass<TowerInfoContainer>(topNode, RawTowerNodeName);
  TowerInfoContainer *_calib_towers = findNode::getClass<TowerInfoContainer>(topNode, CalibTowerNodeName);
  unsigned int ntowers = _raw_towers->size();

  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    unsigned int key = _raw_towers->encode_key(channel);
    TowerInfo *caloinfo_raw = _raw_towers->get_tower_at_channel(channel);
    _calib_towers->get_tower_at_channel(channel)->copy_tower(caloinfo_raw);
    float raw_amplitude = caloinfo_raw->get_energy();
    float calibconst = cdbttree->GetFloatValue(key, m_fieldname);
    _calib_towers->get_tower_at_channel(channel)->set_energy(raw_amplitude * calibconst);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerCalib::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in RawTowerCalibration::CreateNodes");
  }
  // detector node

  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", m_detector));

  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_detector);
    dstNode->addNode(DetNode);
  }

  // towers
  RawTowerNodeName = m_inputNodePrefix + m_detector;
  TowerInfoContainer *_raw_towers = findNode::getClass<TowerInfoContainer>(dstNode, RawTowerNodeName);
  if (!_raw_towers)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + RawTowerNodeName + " node in RawTowerCalibration::CreateNodes");
  }

  CalibTowerNodeName = m_outputNodePrefix + m_detector;
  TowerInfoContainer *_calib_towers = findNode::getClass<TowerInfoContainer>(dstNode, CalibTowerNodeName);
  if (!_calib_towers)
  {
    _calib_towers = dynamic_cast<TowerInfoContainer *>(_raw_towers->CloneMe());
  }
  PHIODataNode<PHObject> *calibtowerNode = new PHIODataNode<PHObject>(_calib_towers, CalibTowerNodeName, "PHObject");
  DetNode->addNode(calibtowerNode);
  return;
}
