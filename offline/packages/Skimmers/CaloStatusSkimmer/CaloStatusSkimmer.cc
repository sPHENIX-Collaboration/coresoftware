#include "CaloStatusSkimmer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// Tower stuff
#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfo.h>
// #include <calobase/TowerInfov3.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// ROOT stuff
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>

// for cluster vertex correction
#include <CLHEP/Vector/ThreeVector.h>
#include <TLorentzVector.h>
#include <array>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

//____________________________________________________________________________..
CaloStatusSkimmer::CaloStatusSkimmer(const std::string &name)
    : SubsysReco(name)
{
  n_eventcounter = 0;
  n_skimcounter = 0;
  n_notowernodecounter = 0;
  std::cout << "CaloStatusSkimmer::CaloStatusSkimmer(const std::string &name) ""Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloStatusSkimmer::~CaloStatusSkimmer()
{
  // std::cout << "CaloStatusSkimmer::~CaloStatusSkimmer() Calling dtor" <<
  // std::endl;
}

//____________________________________________________________________________..
int CaloStatusSkimmer::Init(PHCompositeNode *topNode)
{
  std::cout << "CaloStatusSkimmer::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloStatusSkimmer::process_event(PHCompositeNode *topNode)
{
  n_eventcounter++;
  if (b_do_skim_EMCal)
  {
    TowerInfoContainer *towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
    if (!towers)
    {
      n_notowernodecounter++;
      if (Verbosity > 0)
        std::cout << PHWHERE << "calostatuscheck::process_event: missing TOWERS_CEMC" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const UInt_t ntowers = towers->size();
    uint16_t notinstr_count = 0;
    for (UInt_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_count;
      }
    }
    if (Verbosity > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in EMCal = " << ntowers << ", not-instrumented(empty/missing pckt) towers in EMCal = " << notinstr_count << std::endl;
    }

    if (notinstr_count >= m_EMC_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (b_do_skim_HCal)
  {
    TowerInfoContainer *hcalin_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
    TowerInfoContainer *hcalout_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALOUT");
    if (!hcalin_towers || !hcalout_towers)
    {
      n_notowernodecounter++;
      if (Verbosity > 0)
        std::cout << PHWHERE << "calostatuscheck::process_event: missing TOWERS_HCALIN or TOWERS_HCALOUT" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    const UInt_t ntowers_hcalin = hcalin_towers->size();
    uint16_t notinstr_count_hcalin = 0;
    for (UInt_t ch = 0; ch < ntowers_hcalin; ++ch)
    {
      TowerInfo *tower_in = hcalin_towers->get_tower_at_channel(ch);
      if (tower_in->get_isNotInstr())
      {
        ++notinstr_count_hcalin;
      }
    }

    const UInt_t ntowers_hcalout = hcalout_towers->size();
    uint16_t notinstr_count_hcalout = 0;
    for (UInt_t ch = 0; ch < ntowers_hcalout; ++ch)
    {
      TowerInfo *tower_out = hcalout_towers->get_tower_at_channel(ch);
      if (tower_out->get_isNotInstr())
      {
        ++notinstr_count_hcalout;
      }
    }

    if (Verbosity > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in HCalIn = " << ntowers_hcalin << ", not-instrumented(empty/missing pckt) towers in HCalIn = " << notinstr_count_hcalin << ", ntowers in HCalOut = " << ntowers_hcalout << ", not-instrumented(empty/missing pckt) towers in HCalOut = " << notinstr_count_hcalout << std::endl;
    }

    if (notinstr_count_hcalin >= m_HCal_skim_threshold ||
        notinstr_count_hcalout >= m_HCal_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (b_do_skim_sEPD)
  {
    TowerInfoContainer *sepd_towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_SEPD");
    if (!sepd_towers)
    {
      n_notowernodecounter++;
      if (Verbosity > 0)
        std::cout << PHWHERE << "calostatuscheck::process_event: missing TOWERS_SEPD" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const UInt_t ntowers = sepd_towers->size();
    uint16_t notinstr_count = 0;
    for (UInt_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = sepd_towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_count;
      }
    }

    if (Verbosity > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in sEPD = " << ntowers << ", not-instrumented(empty/missing pckt) towers in sEPD = " << notinstr_count << std::endl;
    }

    if (notinstr_count >= m_sEPD_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (b_do_skim_ZDC)
  {
    TowerInfoContainer *zdc_towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
    if (!zdc_towers)
    {
      n_notowernodecounter++;
      if (Verbosity > 0)
        std::cout << PHWHERE << "calostatuscheck::process_event: missing TOWERS_ZDC" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const UInt_t ntowers = zdc_towers->size();
    uint16_t notinstr_count = 0;
    for (UInt_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = zdc_towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_count;
      }
    }

    if (Verbosity > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in ZDC = " << ntowers << ", not-instrumented(empty/missing pckt) towers in ZDC = " << notinstr_count << std::endl;
    }

    if (notinstr_count >= m_ZDC_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloStatusSkimmer::End(PHCompositeNode *topNode)
{
  std::cout << "CaloStatusSkimmer::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  std::cout << "Total events processed: " << n_eventcounter << std::endl;
  std::cout << "Total events skimmed: " << n_skimcounter << std::endl;
  std::cout << "Total events with missing tower nodes: " << n_notowernodecounter << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
