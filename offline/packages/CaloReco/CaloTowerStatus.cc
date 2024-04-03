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

#include <cmath>
#include <cstdlib>    // for exit
#include <exception>  // for exception
#include <iostream>   // for operator<<, basic_ostream
#include <stdexcept>  // for runtime_error

//____________________________________________________________________________..
CaloTowerStatus::CaloTowerStatus(const std::string &name)
  : SubsysReco(name)
{
  if (Verbosity() > 0)
  {
    std::cout << "CaloTowerStatus::CaloTowerStatus(const std::string &name) Calling ctor" << std::endl;
  }
}

//____________________________________________________________________________..
CaloTowerStatus::~CaloTowerStatus()
{
  delete m_cdbttree_chi2;
  delete m_cdbttree_time;
  if (Verbosity() > 0)
  {
    std::cout << "CaloTowerStatus::~CaloTowerStatus() Calling dtor" << std::endl;
  }
}

//____________________________________________________________________________..
int CaloTowerStatus::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);

  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_detector = "CEMC";
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    m_detector = "HCALIN";
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
  {
    m_detector = "HCALOUT";
  }
  else if (m_dettype == CaloTowerDefs::ZDC)
  {
    m_detector = "ZDC";
  }

  else if (m_dettype == CaloTowerDefs::SEPD)
  {
    m_detector = "SEPD";
  }

  m_calibName_chi2 = m_detector + "_hotTowers_fracBadChi2";
  m_fieldname_chi2 = "fraction";

  std::string calibdir = CDBInterface::instance()->getUrl(m_calibName_chi2);
  if (!calibdir.empty())
  {
    m_cdbttree_chi2 = new CDBTTree(calibdir);
  }
  else
  {
    std::cout << "CaloTowerStatus::::InitRun No masking file for domain " << m_calibName_chi2 << " found, not doing isHot" << std::endl;
    m_doHotChi2 = false;
  }

  m_calibName_time = m_detector + "_meanTime";
  m_fieldname_time = "time";

  calibdir = CDBInterface::instance()->getUrl(m_calibName_time);
  if (!calibdir.empty())
  {
    m_cdbttree_time = new CDBTTree(calibdir);
  }
  else
  {
    std::cout << "CaloTowerStatus::::InitRun no timing info, " << m_calibName_time << " not found, not doing isHot" << std::endl;
    m_doTime = false;
  }

  m_calibName_hotMap = m_detector + "nome";
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_calibName_hotMap = m_detector + "_BadTowerMap";
  }
  m_fieldname_hotMap = "status";

  calibdir = CDBInterface::instance()->getUrl(m_calibName_hotMap);
  if (!calibdir.empty())
  {
    m_cdbttree_hotMap = new CDBTTree(calibdir);
  }
  else
  {
    std::cout << "CaloTowerStatus::::InitRun hot map info, " << m_calibName_hotMap << " not found, not doing isHot" << std::endl;
    m_doHotMap = false;
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
  if (Verbosity() > 0)
  {
    topNode->print();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloTowerStatus::process_event(PHCompositeNode * /*topNode*/)
{
  unsigned int ntowers = m_raw_towers->size();
  float fraction_badChi2 = 0;
  float mean_time = 0;
  int hotMap_val = 0;
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    unsigned int key = m_raw_towers->encode_key(channel);
    m_raw_towers->get_tower_at_channel(channel)->set_status(0);  // resetting status

    if (m_doHotChi2)
    {
      fraction_badChi2 = m_cdbttree_chi2->GetFloatValue(key, m_fieldname_chi2);
    }
    if (m_doTime)
    {
      mean_time = m_cdbttree_time->GetFloatValue(key, m_fieldname_time);
    }
    if (m_doHotMap)
    {
      hotMap_val = m_cdbttree_hotMap->GetIntValue(key, m_fieldname_hotMap);
    }
    float chi2 = m_raw_towers->get_tower_at_channel(channel)->get_chi2();
    float time = m_raw_towers->get_tower_at_channel(channel)->get_time_float();

    if (fraction_badChi2 > fraction_badChi2_threshold && m_doHotChi2)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isHot(true);
    }
    if (std::fabs(time - mean_time) > time_cut && m_doTime)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isBadTime(true);
    }
    if (hotMap_val != 0 && m_doHotMap)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isHot(true);
    }
    if (chi2 > badChi2_treshold)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isBadChi2(true);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerStatus::CreateNodeTree(PHCompositeNode *topNode)
{
  std::string RawTowerNodeName = m_inputNodePrefix + m_detector;
  m_raw_towers = findNode::getClass<TowerInfoContainer>(topNode, RawTowerNodeName);
  if (!m_raw_towers)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + RawTowerNodeName + " node in CaloTowerStatus::CreateNodes");
  }

  return;
}
