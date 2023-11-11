#include "CaloTowerStatus.h"
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
CaloTowerStatus::CaloTowerStatus(const std::string &name)
  : SubsysReco(name)
  , m_dettype(CaloTowerDefs::HCALOUT)
  , m_detector("HCALOUT")
  , m_DETECTOR(TowerInfoContainer::HCAL)
  , m_fieldname("")
  , m_runNumber(-1)
{
  if (Verbosity() > 0) std::cout << "CaloTowerStatus::CaloTowerStatus(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloTowerStatus::~CaloTowerStatus()
{
  delete cdbttree;
  if (Verbosity() > 0) std::cout << "CaloTowerStatus::~CaloTowerStatus() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloTowerStatus::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);

  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_detector = "CEMC";
    std::string default_time_independent_calib = "cemc_pi0_twrSlope_v1_default";

    if (!m_overrideCalibName)
    {
      m_calibName = m_detector + "_BadTowerMap";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "status";
    }
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerStatus::::InitRun No EMCal deadmap found" << std::endl;
      m_doHot = 0;
    }
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    m_detector = "HCALIN";

    if (!m_overrideCalibName)
    {
      m_calibName = m_detector + "_BadTowerMap";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "status";
    }
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerStatus::::InitRun No masking file for domain " << m_calibName << " found, not doing isHot" << std::endl;
      m_doHot = 0;
    }
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
  {
    m_detector = "HCALOUT";

    if (!m_overrideCalibName)
    {
      m_calibName = m_detector + "_BadTowerMap";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "status";
    }
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerStatus::::InitRun No masking file for domain " << m_calibName << " found, not doing isHot" << std::endl;
      m_doHot = 0;
    }
  }
  else if (m_dettype == CaloTowerDefs::ZDC)
  {
    m_detector = "ZDC";

    if (!m_overrideCalibName)
    {
      m_calibName = "noMasker";
    }
    if (!m_overrideFieldName)
    {
      m_fieldname = "noMasker";
    }
    std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
    if (!calibdir.empty())
    {
      cdbttree = new CDBTTree(calibdir);
    }
    else
    {
      std::cout << "CaloTowerStatus::::InitRun No calibration file for domain " << m_calibName << " found" << std::endl;
      exit(1);
    }
  }

  else if (m_dettype == CaloTowerDefs::SEPD)
  {
    m_detector = "SEPD";

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
      std::cout << "CaloTowerStatus::::InitRun No calibration file for domain " << m_calibName << " found" << std::endl;
      exit(1);
    }
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
int CaloTowerStatus::process_event(PHCompositeNode * /*topNode*/)
{
  unsigned int ntowers = _raw_towers->size();
  int mask_status = 0;
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    if (m_doHot) mask_status = cdbttree->GetIntValue(channel, m_fieldname);
    float chi2 = _raw_towers->get_tower_at_channel(channel)->get_chi2();
    if (mask_status > 1 && m_doHot)
    {
      _raw_towers->get_tower_at_channel(channel)->set_isHot(true);
    }
    if (chi2 > 3e3)
    {
      _raw_towers->get_tower_at_channel(channel)->set_isBadChi2(true);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerStatus::CreateNodeTree(PHCompositeNode *topNode)
{
  RawTowerNodeName = m_inputNodePrefix + m_detector;
  _raw_towers = findNode::getClass<TowerInfoContainer>(topNode, RawTowerNodeName);
  if (!_raw_towers)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + RawTowerNodeName + " node in CaloTowerStatus::CreateNodes");
  }

  return;
}
