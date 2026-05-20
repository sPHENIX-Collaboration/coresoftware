#include "CaloTowerStatus.h"
#include "CaloTowerDefs.h"

#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainer.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <iostream>  // for operator<<, basic_ostream

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
int CaloTowerStatus::InitRun(PHCompositeNode *topNode)
{
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

  CreateNodeTree(topNode);

  CDBTTree *cdbttree_chi2{nullptr};

  m_calibName_chi2 = m_detector + "_hotTowers_fracBadChi2";
  m_fieldname_chi2 = "fraction";

  if (!m_directURL_chi2.empty())
  {
    std::cout << "CaloTowerStatus::InitRun: Using direct URL override for chi2: " << m_directURL_chi2 << std::endl;
    cdbttree_chi2 = new CDBTTree(m_directURL_chi2);
  }
  else
  {
    std::string calibdir_chi2 = CDBInterface::instance()->getUrl(m_calibName_chi2);
    if (!calibdir_chi2.empty())
    {
      cdbttree_chi2 = new CDBTTree(calibdir_chi2);
      if (Verbosity() > 0)
      {
        std::cout << "CaloTowerStatus::InitRun Found " << m_calibName_chi2 << "  Doing isHot for frac bad chi2" << std::endl;
      }
    }
    else
    {
      if (m_doAbortNoChi2)
      {
        std::cout << "CaloTowerStatus::InitRun: No chi2 calibration found for " << m_calibName_chi2 << " and abort mode is set. Exiting." << std::endl;
        gSystem->Exit(1);
      }
      m_doHotChi2 = false;
      if (Verbosity() > 0)
      {
        std::cout << "CaloTowerStatus::InitRun No masking file for domain " << m_calibName_chi2 << " found, not doing isHot from isBadChi2" << std::endl;
      }
    }
  }

  CDBTTree *cdbttree_hotMap = nullptr;

  m_calibName_hotMap = m_detector + "_BadTowerMap";
  m_fieldname_hotMap = "status";
  m_fieldname_z_score = m_detector + "_sigma";

  std::string calibdir_hotMap;
  if (!m_directURL_hotMap.empty())
  {
    calibdir_hotMap = m_directURL_hotMap;
    std::cout << "CaloTowerStatus::InitRun: Using direct URL override for hot map: " << calibdir_hotMap << std::endl;
    cdbttree_hotMap = new CDBTTree(calibdir_hotMap);
  }
  else
  {
    calibdir_hotMap = CDBInterface::instance()->getUrl(m_calibName_hotMap);
    if (!calibdir_hotMap.empty())
    {
      cdbttree_hotMap = new CDBTTree(calibdir_hotMap);
      if (Verbosity() > 1)
      {
        std::cout << "CaloTowerStatus::Init " << m_detector << "  hot map found " << m_calibName_hotMap << " Doing isHot" << std::endl;
      }
    }
    else
    {
      if (m_doAbortNoHotMap)
      {
        std::cout << "CaloTowerStatus::InitRun: No hot map found for " << m_calibName_hotMap << " and abort mode is set. Exiting." << std::endl;
        gSystem->Exit(1);
      }
      m_doHotMap = false;
      if (Verbosity() > 1)
      {
        std::cout << "CaloTowerStatus::InitRun hot map info, " << m_calibName_hotMap << " not found, not doing isHot" << std::endl;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTowerStatus::Init " << m_detector << "  doing hotBadChi2=" << std::boolalpha << m_doHotChi2 << " doing hot map=" << std::boolalpha << m_doHotMap << std::endl;
  }

  LoadCalib(cdbttree_chi2, cdbttree_hotMap);

  delete cdbttree_chi2;
  delete cdbttree_hotMap;

  if (Verbosity() > 0)
  {
    topNode->print();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerStatus::LoadCalib(CDBTTree *cdbttree_chi2, CDBTTree *cdbttree_hotMap)
{
  unsigned int ntowers = m_raw_towers->size();
  m_cdbInfo_vec.resize(ntowers);

  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    unsigned int key = m_raw_towers->encode_key(channel);

    if (m_doHotChi2 && cdbttree_chi2)
    {
      m_cdbInfo_vec[channel].fraction_badChi2 = cdbttree_chi2->GetFloatValue(key, m_fieldname_chi2);
    }
    if (m_doHotMap && cdbttree_hotMap)
    {
      m_cdbInfo_vec[channel].hotMap_val = cdbttree_hotMap->GetIntValue(key, m_fieldname_hotMap);
      m_cdbInfo_vec[channel].z_score = cdbttree_hotMap->GetFloatValue(key, m_fieldname_z_score);
    }
  }
}

//____________________________________________________________________________..
int CaloTowerStatus::process_event(PHCompositeNode * /*topNode*/)
{
  unsigned int ntowers = m_raw_towers->size();
  float fraction_badChi2 = 0;
  int hotMap_val = 0;
  float z_score = 0;
  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    // only reset what we will set
    m_raw_towers->get_tower_at_channel(channel)->set_isHot(false);
    m_raw_towers->get_tower_at_channel(channel)->set_isBadChi2(false);

    if (m_doHotChi2)
    {
      fraction_badChi2 = m_cdbInfo_vec[channel].fraction_badChi2;
    }
    if (m_doHotMap)
    {
      hotMap_val = m_cdbInfo_vec[channel].hotMap_val;
      z_score = m_cdbInfo_vec[channel].z_score;
    }
    float chi2 = m_raw_towers->get_tower_at_channel(channel)->get_chi2();
    float adc = m_raw_towers->get_tower_at_channel(channel)->get_energy();

    if (fraction_badChi2 > fraction_badChi2_threshold && m_doHotChi2)
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isHot(true);
    }
    if (m_doHotMap)
    {
      bool is_hot_tower = false;

      // 1. Default behavior: rely on valid positive hotMap status codes only
      if (z_score_threshold == z_score_threshold_default)
      {
        is_hot_tower = (hotMap_val > 0);
      }
      // 2. Custom behavior: evaluate based on the custom z_score threshold
      else
      {
        bool is_dead = (hotMap_val == 1);
        bool exceeds_zscore_limit = (std::abs(z_score) > z_score_threshold);                      // Captures both hot and cold by sigma
        bool is_low_yield_cold = (hotMap_val == 3 && z_score >= -1 * z_score_threshold_default);  // Captures the mean-based cold towers

        is_hot_tower = (is_dead || exceeds_zscore_limit || is_low_yield_cold);
      }

      // Apply the result
      if (is_hot_tower)
      {
        m_raw_towers->get_tower_at_channel(channel)->set_isHot(true);
      }
    }
    if (chi2 > std::min(std::max(badChi2_treshold_const, adc * adc * badChi2_treshold_quadratic), badChi2_treshold_max))
    {
      m_raw_towers->get_tower_at_channel(channel)->set_isBadChi2(true);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerStatus::CreateNodeTree(PHCompositeNode *topNode)
{
  std::string RawTowerNodeName = m_inputNodePrefix + m_detector;
  if (!m_inputNode.empty())
  {
    RawTowerNodeName = m_inputNode;
  }
  m_raw_towers = findNode::getClass<TowerInfoContainer>(topNode, RawTowerNodeName);
  if (!m_raw_towers)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, exiting!"
              << std::endl;
    gSystem->Exit(1);
  }

  return;
}
