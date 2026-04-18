#include "CaloStatusSkimmer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <qautils/QAHistManagerDef.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <TH1.h>

#include <cassert>
#include <cstdint>
#include <iostream>

//____________________________________________________________________________..
CaloStatusSkimmer::CaloStatusSkimmer(const std::string &name)
    : SubsysReco(name)
{
  //std::cout << "CaloStatusSkimmer::CaloStatusSkimmer(const std::string &name) ""Calling ctor" << std::endl;
}


//____________________________________________________________________________..
int CaloStatusSkimmer::Init([[maybe_unused]] PHCompositeNode *topNode)
{
  // std::cout << "CaloStatusSkimmer::Init(PHCompositeNode *topNode) This is Init..." << std::endl;

  if (b_produce_QA_histograms)
  {
    auto* hm = QAHistManagerDef::getHistoManager();
    assert(hm);

    h_EMC_nTowers_notinstr = new TH1F("h_EMC_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in EMCal; nNotInstrTowers; Counts", 24577, -0.5, 24576.5);
    h_EMC_nTowers_notinstr->SetDirectory(nullptr);
    h_HCal_nTowers_notinstr = new TH1F("h_HCal_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in HCal; nNotInstrTowers; Counts", 1537, -0.5, 1536.5);
    h_HCal_nTowers_notinstr->SetDirectory(nullptr);
    h_sEPD_nTowers_notinstr = new TH1F("h_sEPD_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in sEPD; nNotInstrTowers; Counts", 745, -0.5, 744.5);
    h_sEPD_nTowers_notinstr->SetDirectory(nullptr);
    h_ZDC_nTowers_notinstr = new TH1F("h_ZDC_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in ZDC; nNotInstrTowers; Counts", 53, -0.5, 52.5);
    h_ZDC_nTowers_notinstr->SetDirectory(nullptr);

    h_calo_nEvents = new TH1F("h_calo_nEvents", "Number of events", 7, 0.5, 7.5);
    h_calo_nEvents->GetXaxis()->SetBinLabel(1, "Total events processed");
    h_calo_nEvents->GetXaxis()->SetBinLabel(2, "Total events skimmed");
    h_calo_nEvents->GetXaxis()->SetBinLabel(3, "EMCal above not-instr threshold");
    h_calo_nEvents->GetXaxis()->SetBinLabel(4, "HCal above not-instr threshold");
    h_calo_nEvents->GetXaxis()->SetBinLabel(5, "sEPD above not-instr threshold");
    h_calo_nEvents->GetXaxis()->SetBinLabel(6, "ZDC above not-instr threshold");
    h_calo_nEvents->GetXaxis()->SetBinLabel(7, "No TowerInfo nodes found");
    h_calo_nEvents->SetDirectory(nullptr);

    hm->registerHisto(h_calo_nEvents);

    hm->registerHisto(h_EMC_nTowers_notinstr);
    hm->registerHisto(h_HCal_nTowers_notinstr);
    hm->registerHisto(h_sEPD_nTowers_notinstr);
    hm->registerHisto(h_ZDC_nTowers_notinstr);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloStatusSkimmer::process_event(PHCompositeNode *topNode)
{
  n_eventcounter++;
  uint16_t notinstr_EMC = 0;
  uint16_t notinstr_HCalin = 0;
  uint16_t notinstr_HCalout = 0;
  uint16_t notinstr_sEPD = 0;
  uint16_t notinstr_ZDC = 0;

  if (m_EMC_skim_threshold > 0)
  {
    TowerInfoContainer *towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
    if (!towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_CEMC" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const uint32_t ntowers = towers->size();
    for (uint32_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_EMC;
      }
    }
    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in EMCal = " << ntowers << ", not-instrumented(empty/missing pckt) towers in EMCal = " << notinstr_EMC << std::endl;
    }

    if (b_produce_QA_histograms)
    {
      h_EMC_nTowers_notinstr->Fill(notinstr_EMC);
    }

    if (notinstr_EMC >= m_EMC_skim_threshold)
    {
      EMC_skim_count++;
    }
  }

  if (m_HCal_skim_threshold > 0)
  {
    TowerInfoContainer *hcalin_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
    TowerInfoContainer *hcalout_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALOUT");
    if (!hcalin_towers || !hcalout_towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_HCALIN or TOWERS_HCALOUT" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    const uint32_t ntowers_hcalin = hcalin_towers->size();
    for (uint32_t ch = 0; ch < ntowers_hcalin; ++ch)
    {
      TowerInfo *tower_in = hcalin_towers->get_tower_at_channel(ch);
      if (tower_in->get_isNotInstr())
      {
        ++notinstr_HCalin;
      }
    }

    const uint32_t ntowers_hcalout = hcalout_towers->size();
    for (uint32_t ch = 0; ch < ntowers_hcalout; ++ch)
    {
      TowerInfo *tower_out = hcalout_towers->get_tower_at_channel(ch);
      if (tower_out->get_isNotInstr())
      {
        ++notinstr_HCalout;
      }
    }

    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in HCalIn = " << ntowers_hcalin << ", not-instrumented(empty/missing pckt) towers in HCalIn = " << notinstr_HCalin << ", ntowers in HCalOut = " << ntowers_hcalout << ", not-instrumented(empty/missing pckt) towers in HCalOut = " << notinstr_HCalout << std::endl;
    }

    if (b_produce_QA_histograms)
    {
      h_HCal_nTowers_notinstr->Fill(notinstr_HCalin);
      h_HCal_nTowers_notinstr->Fill(notinstr_HCalout);
    }

    if (notinstr_HCalin >= m_HCal_skim_threshold ||
        notinstr_HCalout >= m_HCal_skim_threshold)
    {
      HCal_skim_count++;
    }
  }

  if (m_sEPD_skim_threshold > 0)
  {
    TowerInfoContainer *sepd_towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_SEPD");
    if (!sepd_towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_SEPD" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const uint32_t ntowers = sepd_towers->size();
    for (uint32_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = sepd_towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_sEPD;
      }
    }

    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in sEPD = " << ntowers << ", not-instrumented(empty/missing pckt) towers in sEPD = " << notinstr_sEPD << std::endl;
    }

    if(b_produce_QA_histograms)
    {
      h_sEPD_nTowers_notinstr->Fill(notinstr_sEPD);
    }

    if (notinstr_sEPD >= m_sEPD_skim_threshold)
    {
      sEPD_skim_count++;
    }
  }

  if (m_ZDC_skim_threshold > 0)
  {
    TowerInfoContainer *zdc_towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
    if (!zdc_towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_ZDC" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const uint32_t ntowers = zdc_towers->size();
    for (uint32_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = zdc_towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_ZDC;
      }
    }

    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in ZDC = " << ntowers << ", not-instrumented(empty/missing pckt) towers in ZDC = " << notinstr_ZDC << std::endl;
    }

    if(b_produce_QA_histograms)
    {
      h_ZDC_nTowers_notinstr->Fill(notinstr_ZDC);
    }

    if (notinstr_ZDC >= m_ZDC_skim_threshold)
    {
      ZDC_skim_count++;
    }
  }

  // If any of the enabled skimming conditions are met, then increment the skim counter and return ABORTEVENT to skip the event
  if ((m_EMC_skim_threshold > 0 && notinstr_EMC >= m_EMC_skim_threshold) || (m_HCal_skim_threshold > 0 && (notinstr_HCalin >= m_HCal_skim_threshold || notinstr_HCalout >= m_HCal_skim_threshold)) || (m_sEPD_skim_threshold > 0 && notinstr_sEPD >= m_sEPD_skim_threshold) || (m_ZDC_skim_threshold > 0 && notinstr_ZDC >= m_ZDC_skim_threshold))
  {
    n_skimcounter++;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloStatusSkimmer::End([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "CaloStatusSkimmer::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  std::cout << "CaloStatusSkimmer::End Total events processed: " << n_eventcounter << std::endl;
  std::cout << "CaloStatusSkimmer::End Total events skimmed: " << n_skimcounter << std::endl;
  std::cout << "CaloStatusSkimmer::End Total events with missing tower nodes: " << n_notowernodecounter << std::endl;

  if (b_produce_QA_histograms)
  {
    h_calo_nEvents->SetBinContent(1, n_eventcounter);
    h_calo_nEvents->SetBinContent(2, n_skimcounter);
    h_calo_nEvents->SetBinContent(3, EMC_skim_count);
    h_calo_nEvents->SetBinContent(4, HCal_skim_count);
    h_calo_nEvents->SetBinContent(5, sEPD_skim_count);
    h_calo_nEvents->SetBinContent(6, ZDC_skim_count);
    h_calo_nEvents->SetBinContent(7, n_notowernodecounter);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
