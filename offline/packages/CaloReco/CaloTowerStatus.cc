#include "CaloTowerStatus.h"
#include "CaloTowerDefs.h"

#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoDefs.h>
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
  if (Verbosity() > 0)
  {
    std::cout << "CaloTowerStatus::~CaloTowerStatus() Calling dtor" << std::endl;
  }
  delete m_cdbttree_chi2;
  delete m_cdbttree_time;
  delete m_cdbttree_hotMap;
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
    if (Verbosity() > 0)
    {
       std::cout << "CaloTowerStatus::InitRun Found " << m_calibName_chi2 << "  Doing isHot for frac bad chi2" << std::endl;
    }
  }
  else
  {
    if (use_directURL_chi2)
    {
      calibdir = m_directURL_chi2;
      std::cout << "CaloTowerStatus::InitRun: Using default hotBadChi2" << calibdir << std::endl;
      m_cdbttree_chi2 = new CDBTTree(calibdir);
    }
    else 
    {
      m_doHotChi2 = false;
      if (Verbosity() > 0)
      {
        std::cout << "CaloTowerStatus::InitRun No masking file for domain " << m_calibName_chi2 << " found, not doing isHot from isBadChi2" << std::endl;
      }
    }
  }

  m_calibName_time = m_detector + "_meanTime";
  m_fieldname_time = "time";

  calibdir = CDBInterface::instance()->getUrl(m_calibName_time);
  if (!calibdir.empty())
  {
    m_cdbttree_time = new CDBTTree(calibdir);
    if (Verbosity() > 0)
    {
      std::cout << "CaloTowerStatus::InitRun Found " << m_calibName_time << " not Doing isBadTime" << std::endl;
    }
  }
  else
  {
    if (use_directURL_time)
    {
      calibdir = m_directURL_time;
      std::cout << "CaloTowerStatus::InitRun: Using default time  " << calibdir << std::endl;
      m_cdbttree_time = new CDBTTree(calibdir);
    }
    else
    {
      m_doTime = false;
      if (Verbosity() > 1)
      {
        std::cout << "CaloTowerStatus::InitRun no timing info, " << m_calibName_time << " not found, not doing isBadTime" << std::endl;
      }
    }
  }

  m_calibName_hotMap = m_detector + "nome";
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_calibName_hotMap = m_detector + "_BadTowerMap";
  }
  m_fieldname_hotMap = "status";
  m_fieldname_z_score = m_detector + "_sigma";

  calibdir = CDBInterface::instance()->getUrl(m_calibName_hotMap);
  if (!calibdir.empty())
  {
    m_cdbttree_hotMap = new CDBTTree(calibdir);
    if (Verbosity() > 1)
    {
      std::cout << "CaloTowerStatus::Init " << m_detector << "  hot map found " << m_calibName_hotMap << " Ddoing isHot" << std::endl;
    }
  }
  else
  {
    if (m_doAbortNoHotMap)
    {
      std::cout << "CaloTowerStatus::InitRun: No hot map.. exiting" << std::endl;
      gSystem->Exit(1);
    }
    if (use_directURL_hotMap)
    {
      calibdir = m_directURL_hotMap;
      std::cout << "CaloTowerStatus::InitRun: Using default map " << calibdir << std::endl;
      m_cdbttree_hotMap = new CDBTTree(calibdir);
    }
    else
    {
      m_doHotMap = false;
      if (Verbosity() > 1)
      {
        std::cout << "CaloTowerStatus::InitRun hot map info, " << m_calibName_hotMap << " not found, not doing isHot" << std::endl;
      }
    }  
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTowerStatus::Init " << m_detector << "  doing time status =" <<  std::boolalpha << m_doTime << "  doing hotBadChi2=" <<  std::boolalpha << m_doHotChi2 << " doing hot map=" << std::boolalpha << m_doHotMap << std::endl;
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

void CaloTowerStatus::emcal_propogate_isBadChi2(const std::vector<std::vector<int>> &badChi2_IB_vec, const std::vector<std::vector<int>> &towers_IB_vec)
{
  unsigned int ntowers = m_raw_towers->size();
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    std::pair<int, int> sector_ib = TowerInfoDefs::getEMCalSectorIB(channel);
    int sector = sector_ib.first;
    int ib = sector_ib.second;
    int badChi2_towers = badChi2_IB_vec[sector][ib];
    int towers = towers_IB_vec[sector][ib];
    // if entire IB is bad then skip
    if (towers == 0)
    {
      continue;
    }
    float badChi2_IB_frac = badChi2_towers * 1. / towers;
    if(badChi2_IB_frac > m_badChi2_IB_threshold)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isBadChi2(true);
    }
  }
}

//____________________________________________________________________________..
int CaloTowerStatus::process_event(PHCompositeNode * /*topNode*/)
{
  unsigned int ntowers = m_raw_towers->size();
  float fraction_badChi2 = 0;
  float mean_time = 0;
  int hotMap_val = 0;
  float z_score = 0;

  std::vector<std::vector<int>> badChi2_IB_vec(emcal_sector, std::vector<int>(emcal_ib_per_sector, 0));
  std::vector<std::vector<int>> towers_IB_vec(emcal_sector, std::vector<int>(emcal_ib_per_sector, 0));

  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    unsigned int key = m_raw_towers->encode_key(channel);
    int sector = -1;
    int ib = -1;

    if(m_dettype == CaloTowerDefs::CEMC)
    {
      std::pair<int, int> sector_ib = TowerInfoDefs::getEMCalSectorIB(channel);
      sector = sector_ib.first;
      ib = sector_ib.second;
    }
    // only reset what we will set
    m_raw_towers->get_tower_at_channel(channel)->set_isHot(false);
    m_raw_towers->get_tower_at_channel(channel)->set_isBadTime(false);
    m_raw_towers->get_tower_at_channel(channel)->set_isBadChi2(false);

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
      if(!m_isSim)
      {
        z_score = m_cdbttree_hotMap->GetFloatValue(key, m_fieldname_z_score);
      }
    }
    float chi2 = m_raw_towers->get_tower_at_channel(channel)->get_chi2();
    float time = m_raw_towers->get_tower_at_channel(channel)->get_time_float();
    float adc = m_raw_towers->get_tower_at_channel(channel)->get_energy();

    if (fraction_badChi2 > fraction_badChi2_threshold && m_doHotChi2)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isHot(true);
    }
    if (!m_raw_towers->get_tower_at_channel(channel)->get_isZS() && std::fabs(time - mean_time) > time_cut && m_doTime)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isBadTime(true);
    }

    bool is_bad_tower_data = hotMap_val == 1 ||                                              // dead
                             std::fabs(z_score) > z_score_threshold ||                       // hot or cold
                             (hotMap_val == 3 && z_score >= -1 * z_score_threshold_default); // cold part 2

    bool is_bad_tower_sim = hotMap_val != 0;

    bool is_bad_tower = (!m_isSim && is_bad_tower_data) || (m_isSim && is_bad_tower_sim);

    if (is_bad_tower && m_doHotMap)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isHot(true);
    }

    if(m_dettype == CaloTowerDefs::CEMC && m_raw_towers->get_tower_at_channel(channel)->get_isGood())
    {
      ++towers_IB_vec[sector][ib];
    }

    if (chi2 > std::min(std::max(badChi2_treshold_const, adc * adc * badChi2_treshold_quadratic),badChi2_treshold_max))
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isBadChi2(true);

      if(m_dettype == CaloTowerDefs::CEMC)
      {
        ++badChi2_IB_vec[sector][ib];
      }
    }
  }

  // propagate the isBadChi2 status to entire interface board if threshold is exceeded
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    emcal_propogate_isBadChi2(badChi2_IB_vec, towers_IB_vec);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerStatus::CreateNodeTree(PHCompositeNode *topNode)
{
  std::string RawTowerNodeName = m_inputNodePrefix + m_detector;
  m_raw_towers = findNode::getClass<TowerInfoContainer>(topNode, RawTowerNodeName);
  if (!m_raw_towers)
  {
    std::cout << Name() << "::" << m_detector.c_str() << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + RawTowerNodeName + " node in CaloTowerStatus::CreateNodes");
  }

  return;
}
