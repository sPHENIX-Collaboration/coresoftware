#include "CaloValid.h"

// Calo includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <qautils/QAHistManagerDef.h>

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/RunHeader.h>  // for runnumber

#include <ffarawobjects/Gl1Packet.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/RunnumberRange.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TSystem.h>

#include <cassert>
#include <cmath>   // log10, pow, sqrt, std::abs, M_PI
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

namespace
{
  // helper: faster alpha calc without pow/sqrt
  inline float clusterAlpha(const float e1, const float e2)
  {
    const float denom = e1 + e2;
    if (denom == 0.F) { return std::numeric_limits<float>::infinity(); }
    return std::abs(e1 - e2) / denom;
  }

  inline bool inTrigOfInterest(const std::vector<int>& trigs, int bit)
  {
    return std::find(trigs.begin(), trigs.end(), bit) != trigs.end();
  }
}  // namespace

CaloValid::CaloValid(const std::string& name)
  : SubsysReco(name)
{
}

CaloValid::~CaloValid()
{
  for (int i = 0; i < 128 * 192; i++)
  {
    delete h_cemc_channel_pedestal[i];
    delete h_cemc_channel_energy[i];
  }
  for (int i = 0; i < 32 * 48; i++)
  {
    delete h_ihcal_channel_pedestal[i];
    delete h_ihcal_channel_energy[i];
    delete h_ohcal_channel_pedestal[i];
    delete h_ohcal_channel_energy[i];
  }
  delete trigAna;
}

int CaloValid::Init(PHCompositeNode* /*unused*/)
{
  if (m_debug)
  {
    std::cout << "In CaloValid::Init" << std::endl;
  }

  createHistos();
  trigAna = new TriggerAnalyzer();
  if (m_debug)
  {
    std::cout << "Leaving CaloValid::Init" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::InitRun(PHCompositeNode* topNode)
{
  RunHeader* runhdr = findNode::getClass<RunHeader>(topNode, "RunHeader");

  if (runhdr)
  {
    const int runnumber = runhdr->get_RunNumber();

    if (runnumber >= RunnumberRange::RUN2PP_FIRST && runnumber <= RunnumberRange::RUN2PP_LAST)
    {
      m_species = "pp";
      if (Verbosity() > 0)
      {
        std::cout << "This run is from Run-2 p+p.\n";
      }
    }
    else if (runnumber >= RunnumberRange::RUN2AUAU_FIRST && runnumber <= RunnumberRange::RUN2AUAU_LAST)
    {
      m_species = "AuAu";
      if (Verbosity() > 0)
      {
        std::cout << "This run is from Run-2 Au+Au.\n";
      }
    }
    else if (runnumber >= RunnumberRange::RUN3AUAU_FIRST && runnumber <= RunnumberRange::RUN3AUAU_LAST)
    {
      m_species = "AuAu";
      if (Verbosity() > 0)
      {
        std::cout << "This run is from Run-3 Au+Au.\n";
      }
    }
    else if (Verbosity() > 0)
    {
      std::cout << "Run number is out of range. Check RunnumberRange.h . Using pp as default.\n";
    }
  }
  else if (Verbosity() > 0)
  {
    std::cout << "No RunHeader node found. Using pp as default species.\n";
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;
  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::process_towers(PHCompositeNode* topNode)
{
  // --------------------------- Event header -------------------------------- //
  EventHeader* eventheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  int event_number = 0;
  if (eventheader)
  {
    if (eventheader->isValid())
    {
      event_number = eventheader->get_EvtSequence();
    }
  }
  else
  {
    std::cout << "GlobalQA::process_event()  No event header" << std::endl;
  }

  if (m_debug)
  {
    std::cout << _eventcounter << std::endl;
  }

  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  float totalcemc = 0.F;
  float totalihcal = 0.F;
  float totalohcal = 0.F;
  float totalmbd = 0.F;

  float emcaldownscale = 0.F;
  float ihcaldownscale = 0.F;
  float ohcaldownscale = 0.F;
  float mbddownscale = 0.F;
  float adc_threshold = 0.F;
  float emcal_hit_threshold = 0.F;
  float emcal_highhit_threshold = 0.F;
  float ohcal_hit_threshold = 0.F;
  float ohcal_highhit_threshold = 0.F;
  float ihcal_hit_threshold = 0.F;
  float ihcal_highhit_threshold = 0.F;

  if (m_species == "AuAu")
  {
    emcaldownscale = 1350000.F / 800.F;
    ihcaldownscale = 55000.F / 300.F;
    ohcaldownscale = 265000.F / 600.F;
    mbddownscale = 2800.F;
    adc_threshold = 15.F;

    emcal_hit_threshold = 0.5F;  // GeV
    ohcal_hit_threshold = 0.5F;
    ihcal_hit_threshold = 0.25F;

    emcal_highhit_threshold = 3.0F;
    ohcal_highhit_threshold = 3.0F;
    ihcal_highhit_threshold = 3.0F;
  }
  else
  {
    emcaldownscale = 100000.F / 800.F;
    ihcaldownscale = 4000.F / 300.F;
    ohcaldownscale = 25000.F / 600.F;
    mbddownscale = 200.F;
    adc_threshold = 100.F;

    emcal_hit_threshold = 0.5F;  // GeV
    ohcal_hit_threshold = 0.5F;
    ihcal_hit_threshold = 0.25F;

    emcal_highhit_threshold = 3.0F;
    ohcal_highhit_threshold = 3.0F;
    ihcal_highhit_threshold = 3.0F;
  }

  // ---------------------------------- vertex -------------------------------- //
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << "CaloValid GlobalVertexMap node is missing" << std::endl;
  }
  float vtx_z = std::numeric_limits<float>::quiet_NaN();

  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex* vtx = vertexmap->begin()->second;
    if (vtx)
    {
      vtx_z = static_cast<float>(vtx->get_z());
    }
    h_vtx_z_raw->Fill(vtx_z);
    if (std::fabs(vtx_z) < 20.0F)
    {
      h_vtx_z_cut->Fill(vtx_z);
    }
  }

  // --------------------------- trigger and GL1 ------------------------------ //
  if (trigAna)
  {
    trigAna->decodeTriggers(topNode);
  }
  else if (m_debug)
  {
    std::cout << "[ERROR] No TriggerAnalyzer pointer!\n";
  }

  std::vector<int> scaledActiveBits;
  scaledActiveBits.reserve(triggerIndices.size());

  for (int bit : triggerIndices)
  {
    if (trigAna->didTriggerFire(bit))
    {
      scaledActiveBits.push_back(bit);
    }
  }

  bool scaledBits[64] = {false};
  uint64_t raw[64] = {0};
  uint64_t live[64] = {0};

  Gl1Packet* gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, 14001);
  if (!gl1PacketInfo)
  {
    gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1PacketInfo)
    {
      std::cout << PHWHERE << "GlobalQA::process_event: GL1Packet node is missing" << std::endl;
    }
  }

  if (gl1PacketInfo)
  {
    uint64_t triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++)
    {
      const bool trig_decision = ((triggervec & 0x1U) == 0x1U);
      scaledBits[i] = trig_decision;

      raw[i] = gl1PacketInfo->lValue(i, 0);
      live[i] = gl1PacketInfo->lValue(i, 1);

      if (trig_decision)
      {
        h_triggerVec->Fill(i);
      }
      triggervec = (triggervec >> 1U) & 0xffffffffU;
    }
  }

  // --------------------------- calibrated towers --------------------------- //
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if (towers)
    {
      const int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        const float offlineenergy = tower->get_energy();
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        h_cemc_e_chi2->Fill(offlineenergy, tower->get_chi2());
        const float timef = tower->get_time();
        h_emcaltime_cut->Fill(timef);
        const bool isGood = tower->get_isGood();
        uint8_t status = tower->get_status();
        h_emcal_tower_e->Fill(offlineenergy);
        h_emcal_tower_e_wide_range->Fill(offlineenergy);
        if (tower->get_isSaturated())
        {
          h_emcal_tower_e_saturated->Fill(offlineenergy);
        }
        const float pedestal = tower->get_pedestal();
        h_cemc_channel_pedestal[channel]->Fill(pedestal);

        for (int is = 0; is < 8; is++)
        {
          if ((status & 1U) != 0U)
          {
            h_cemc_status->Fill(is);
          }
          status = static_cast<uint8_t>(status >> 1U);
        }

        totalcemc += offlineenergy;
        h_emcaltime->Fill(timef);
        if (offlineenergy > emcal_hit_threshold)
        {
          h_cemc_etaphi_time->Fill(ieta, iphi, timef);
          h_cemc_etaphi->Fill(ieta, iphi);
          if (isGood && (scaledBits[10] || scaledBits[11]))
          {
            h_cemc_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
          }
          h_cemc_etaphi_badChi2->Fill(ieta, iphi, tower->get_isBadChi2() ? 1.0 : 0.0);
        }

        if (offlineenergy > 0.25F)
        {
          h_cemc_etaphi_fracHit->Fill(ieta, iphi, 1.0);
        }
        else
        {
          h_cemc_etaphi_fracHit->Fill(ieta, iphi, 0.0);
        }
        if (offlineenergy > emcal_highhit_threshold)
        {
          h_cemc_etaphi_time_highhit->Fill(ieta, iphi, timef);
          h_cemc_etaphi_highhit->Fill(ieta, iphi);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if (towers)
    {
      const int size = towers->size();
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        const float offlineenergy = tower->get_energy();
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);

        const float timef = tower->get_time();
        h_ihcaltime_cut->Fill(timef);
        h_ihcal_e_chi2->Fill(offlineenergy, tower->get_chi2());
        const bool isGood = tower->get_isGood();
        h_ihcal_status->Fill(tower->get_status());
        const float pedestal = tower->get_pedestal();
        h_ihcal_channel_pedestal[channel]->Fill(pedestal);
        h_ihcal_tower_e->Fill(offlineenergy);
        h_ihcal_tower_e_wide_range->Fill(offlineenergy);
        if (tower->get_isSaturated())
        {
          h_ihcal_tower_e_saturated->Fill(offlineenergy);
        }

        uint8_t status = tower->get_status();
        for (int is = 0; is < 8; is++)
        {
          if ((status & 1U) != 0U)
          {
            h_ihcal_status->Fill(is);
          }
          status = static_cast<uint8_t>(status >> 1U);
        }

        totalihcal += offlineenergy;
        h_ihcaltime->Fill(timef);

        if (offlineenergy > ihcal_hit_threshold)
        {
          h_ihcal_etaphi->Fill(ieta, iphi);
          h_ihcal_etaphi_time->Fill(ieta, iphi, timef);
          if (isGood && (scaledBits[10] || scaledBits[11]))
          {
            h_ihcal_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
          }
          h_ihcal_etaphi_badChi2->Fill(ieta, iphi, tower->get_isBadChi2() ? 1.0 : 0.0);
        }
        if (offlineenergy > ihcal_highhit_threshold)
        {
          h_ihcal_etaphi_time_highhit->Fill(ieta, iphi, timef);
          h_ihcal_etaphi_highhit->Fill(ieta, iphi);
        }
      }
    }
  }

  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if (towers)
    {
      const int size = towers->size();
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        const float offlineenergy = tower->get_energy();
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        const float timef = tower->get_time();
        h_ohcaltime_cut->Fill(timef);
        h_ohcal_e_chi2->Fill(offlineenergy, tower->get_chi2());
        const bool isGood = tower->get_isGood();
        h_ohcal_status->Fill(tower->get_status());
        const float pedestal = tower->get_pedestal();
        h_ohcal_channel_pedestal[channel]->Fill(pedestal);
        h_ohcal_tower_e->Fill(offlineenergy);
        h_ohcal_tower_e_wide_range->Fill(offlineenergy);
        if (tower->get_isSaturated())
        {
          h_ohcal_tower_e_saturated->Fill(offlineenergy);
        }

        uint8_t status = tower->get_status();
        for (int is = 0; is < 8; is++)
        {
          if ((status & 1U) != 0U)
          {
            h_ohcal_status->Fill(is);
          }
          status = static_cast<uint8_t>(status >> 1U);
        }

        totalohcal += offlineenergy;
        h_ohcaltime->Fill(timef);

        if (offlineenergy > ohcal_hit_threshold)
        {
          h_ohcal_etaphi->Fill(ieta, iphi);
          h_ohcal_etaphi_time->Fill(ieta, iphi, timef);
          if (isGood && (scaledBits[10] || scaledBits[11]))
          {
            h_ohcal_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
          }
          h_ohcal_etaphi_badChi2->Fill(ieta, iphi, tower->get_isBadChi2() ? 1.0 : 0.0);
        }
        if (offlineenergy > ohcal_highhit_threshold)
        {
          h_ohcal_etaphi_time_highhit->Fill(ieta, iphi, timef);
          h_ohcal_etaphi_highhit->Fill(ieta, iphi);
        }
      }
    }
  }

  // -------------------------- raw tower ------------------------------ //
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
    if (towers)
    {
      const int size = towers->size();
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        const float raw_time = tower->get_time();
        if (tower->get_isZS())
        {
          h_cemc_channel_energy[channel]->Fill(tower->get_energy());
        }

        const float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_cemc_etaphi_fracHitADC->Fill(ieta, iphi, 1.0);
          h_cemc_etaphi_time_raw->Fill(ieta, iphi, raw_time);
        }
        else
        {
          h_cemc_etaphi_fracHitADC->Fill(ieta, iphi, 0.0);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALOUT");
    if (towers)
    {
      const int size = towers->size();
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        const float raw_time = tower->get_time();
        if (tower->get_isZS())
        {
          h_ohcal_channel_energy[channel]->Fill(tower->get_energy());
        }

        const float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_ohcal_etaphi_time_raw->Fill(ieta, iphi, raw_time);
          h_ohcal_etaphi_fracHitADC->Fill(ieta, iphi, 1.0);
        }
        else
        {
          h_ohcal_etaphi_fracHitADC->Fill(ieta, iphi, 0.0);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
    if (towers)
    {
      const int size = towers->size();
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        const float raw_time = tower->get_time();
        if (tower->get_isZS())
        {
          h_ihcal_channel_energy[channel]->Fill(tower->get_energy());
        }

        const float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_ihcal_etaphi_time_raw->Fill(ieta, iphi, raw_time);
          h_ihcal_etaphi_fracHitADC->Fill(ieta, iphi, 1.0);
        }
        else
        {
          h_ihcal_etaphi_fracHitADC->Fill(ieta, iphi, 0.0);
        }
      }
    }
  }

  // --------------------------- MBD ---------------------------------------- //
  MbdPmtContainer* bbcpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!bbcpmts)
  {
    std::cout << "CaloValid::process_event: Could not find MbdPmtContainer," << std::endl;
  }

  int hits = 0;
  if (bbcpmts)
  {
    const int nPMTs = bbcpmts->get_npmt();
    for (int i = 0; i < nPMTs; i++)
    {
      MbdPmtHit* mbdpmt = bbcpmts->get_pmt(i);
      const float pmtadc = mbdpmt->get_q();
      totalmbd += pmtadc;
      if (pmtadc > 0.4F)
      {
        hits++;
      }
    }
  }
  h_mbd_hits->Fill(hits);

  h_emcal_mbd_correlation->Fill(totalcemc / emcaldownscale, totalmbd / mbddownscale);
  h_ihcal_mbd_correlation->Fill(totalihcal / ihcaldownscale, totalmbd / mbddownscale);
  h_ohcal_mbd_correlation->Fill(totalohcal / ohcaldownscale, totalmbd / mbddownscale);
  h_emcal_hcal_correlation->Fill(totalcemc / emcaldownscale, totalohcal / ohcaldownscale);

  // ------------------------------ clusters & pi0 -------------------------- //
  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if (!clusterContainer)
  {
    std::cout << PHWHERE
              << "CaloValid::process_event - Fatal Error - CLUSTERINFO_CEMC node is missing."
              << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // geometry for hot tower/cluster masking (still required by codebase)
  const std::string towergeomnodename = "TOWERGEOM_CEMC";
  RawTowerGeomContainer* m_geometry = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!m_geometry)
  {
    std::cout << Name() << "::CreateNodeTree"
              << ": Could not find node " << towergeomnodename << std::endl;
    gSystem->Exit(1);
  }

  // pT cuts instead of energy cuts
  constexpr float minClusPt1 = 1.5F;
  constexpr float minClusPt2 = 1.0F;
  constexpr float emcMaxClusE = 100.F;
  constexpr float maxAlpha = 0.6F;

  // low activity gate, using 0.1 scale-down factor
  if (totalcemc < 0.1F * emcaldownscale)
  {
    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();

    for (auto clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; ++clusterIter)
    {
      RawCluster* recoCluster = clusterIter->second;
      CLHEP::Hep3Vector vertex(0, 0, 0);
      const CLHEP::Hep3Vector e_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

      const float clusE = static_cast<float>(e_vec_cluster.mag());
      const float clus_eta = static_cast<float>(e_vec_cluster.pseudoRapidity());
      const float clus_phi = static_cast<float>(e_vec_cluster.phi());
      const float clus_pt = static_cast<float>(e_vec_cluster.perp());
      const float clus_chisq = static_cast<float>(recoCluster->get_chi2());

      h_clusE->Fill(clusE);

      // leading photon candidate cut
      if (clus_pt < minClusPt1 || clusE > emcMaxClusE)
      {
        continue;
      }
      if (clus_chisq > 4.F)
      {
        continue;
      }

      h_etaphi_clus->Fill(clus_eta, clus_phi);

      TLorentzVector photon1;
      photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);

      for (auto clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; ++clusterIter2)
      {
        if (clusterIter == clusterIter2)
        {
          continue;
        }

        RawCluster* recoCluster2 = clusterIter2->second;

        CLHEP::Hep3Vector vertex2(0, 0, 0);
        const CLHEP::Hep3Vector e_vec_cluster2 = RawClusterUtility::GetECoreVec(*recoCluster2, vertex2);

        const float clus2E = static_cast<float>(e_vec_cluster2.mag());
        const float clus2_eta = static_cast<float>(e_vec_cluster2.pseudoRapidity());
        const float clus2_phi = static_cast<float>(e_vec_cluster2.phi());
        const float clus2_pt = static_cast<float>(e_vec_cluster2.perp());
        const float clus2_chisq = static_cast<float>(recoCluster2->get_chi2());

        // partner photon candidate cut
        if (clus2_pt < minClusPt2 || clus2E > emcMaxClusE)
        {
          continue;
        }
        if (clus2_chisq > 4.F)
        {
          continue;
        }

        // energy asymmetry cut (alpha)
        if (clusterAlpha(clusE, clus2E) > maxAlpha)
        {
          continue;
        }

        TLorentzVector photon2;
        photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);

        const TLorentzVector pi0 = photon1 + photon2;
        const float pi0Mass = static_cast<float>(pi0.M());
        const unsigned int lt_eta = recoCluster->get_lead_tower().first;
        const unsigned int lt_phi = recoCluster->get_lead_tower().second;

        // integer index for the 8x8 "IB" macro-cell
        const int IB_num = static_cast<int>((lt_eta / 8U) * 32U + (lt_phi / 8U));

        // require vertex |z| < 20 cm for filling mass
        if (std::fabs(vtx_z) < 20.0F)
        {
          for (int bit : scaledActiveBits)
          {
            if (!inTrigOfInterest(triggerIndices, bit))
            {
              continue;
            }
            h_pi0_trigIB_mass->Fill(
                static_cast<double>(bit),
                static_cast<double>(IB_num),
                static_cast<double>(pi0Mass));
          }
          h_InvMass->Fill(pi0Mass);
        }
      }
    }  // end cluster loop
  }

  // ----------------- Trigger / alignment ---------------------------- //
  float leading_cluster_ecore = 0.F;
  float leading_cluster_eta = 0.F;
  float leading_cluster_phi = 0.F;
  const int evtNum_overK = event_number / 1000;

  if (clusterContainer)
  {
    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();

    for (auto clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; ++clusterIter)
    {
      RawCluster* recoCluster = clusterIter->second;
      if (recoCluster->get_chi2() > 2.F)
      {
        continue;
      }

      CLHEP::Hep3Vector vertex(0, 0, 0);
      const CLHEP::Hep3Vector e_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

      const float clusE = static_cast<float>(e_vec_cluster.mag());
      const float clusEta = static_cast<float>(e_vec_cluster.pseudoRapidity());
      const float clusPhi = static_cast<float>(e_vec_cluster.phi());
      if (clusE > leading_cluster_ecore)
      {
        leading_cluster_ecore = clusE;
        leading_cluster_eta = clusEta;
        leading_cluster_phi = clusPhi;
      }
    }

    for (int i = 0; i < 64; i++)
    {
      if (!scaledBits[i])
      {
        continue;
      }

      pr_ldClus_trig->Fill(i, leading_cluster_ecore);

      if (!inTrigOfInterest(trigOfInterest, i))
      {
        continue;
      }

      h_edist[i]->Fill(leading_cluster_eta, leading_cluster_phi);
      h_ldClus_trig[i]->Fill(leading_cluster_ecore);
      pr_evtNum_ldClus_trig[i]->Fill(evtNum_overK, leading_cluster_ecore);
      if (raw[i] > 0U)
      {
        pr_rejection[i]->Fill(evtNum_overK,
                              static_cast<float>(raw[10]) / static_cast<float>(raw[i]));
        pr_livetime[i]->Fill(evtNum_overK,
                             static_cast<float>(live[i]) / static_cast<float>(raw[i]));
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::End(PHCompositeNode* topNode)
{
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  //------EmCal-----//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if (towers)
    {
      const int size = towers->size();

      auto* h_CaloValid_cemc_etaphi_pedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(std::format("{}cemc_etaphi_pedRMS", getHistoPrefix())));
      auto* h_CaloValid_cemc_etaphi_ZSpedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(std::format("{}cemc_etaphi_ZSpedRMS", getHistoPrefix())));

      for (int channel = 0; channel < size; channel++)
      {
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        const float ped_rms = static_cast<float>(h_cemc_channel_pedestal[channel]->GetRMS());
        h_CaloValid_cemc_etaphi_pedRMS->Fill(ieta, iphi, ped_rms);
        MirrorHistogram(h_cemc_channel_energy[channel]);
        const double rmsZS = h_cemc_channel_energy[channel]->GetRMS();
        h_CaloValid_cemc_etaphi_ZSpedRMS->Fill(ieta, iphi, rmsZS);
      }
    }
  }
  //------IHCal------//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if (towers)
    {
      const int size = towers->size();

      auto* h_CaloValid_ihcal_etaphi_pedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(std::format("{}ihcal_etaphi_pedRMS", getHistoPrefix())));
      auto* h_CaloValid_ihcal_etaphi_ZSpedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(std::format("{}ihcal_etaphi_ZSpedRMS", getHistoPrefix())));

      for (int channel = 0; channel < size; channel++)
      {
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        const float ped_rms = static_cast<float>(h_ihcal_channel_pedestal[channel]->GetRMS());
        h_CaloValid_ihcal_etaphi_pedRMS->Fill(ieta, iphi, ped_rms);
        MirrorHistogram(h_ihcal_channel_energy[channel]);
        const double rmsZS = h_ihcal_channel_energy[channel]->GetRMS();
        h_CaloValid_ihcal_etaphi_ZSpedRMS->Fill(ieta, iphi, rmsZS);
      }
    }
  }
  //------OHCal-----//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if (towers)
    {
      const int size = towers->size();

      auto* h_CaloValid_ohcal_etaphi_pedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(std::format("{}ohcal_etaphi_pedRMS", getHistoPrefix())));
      auto* h_CaloValid_ohcal_etaphi_ZSpedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(std::format("{}ohcal_etaphi_ZSpedRMS", getHistoPrefix())));

      for (int channel = 0; channel < size; channel++)
      {
        const unsigned int towerkey = towers->encode_key(channel);
        const int ieta = towers->getTowerEtaBin(towerkey);
        const int iphi = towers->getTowerPhiBin(towerkey);
        const float ped_rms = static_cast<float>(h_ohcal_channel_pedestal[channel]->GetRMS());
        h_CaloValid_ohcal_etaphi_pedRMS->Fill(ieta, iphi, ped_rms);
        MirrorHistogram(h_ohcal_channel_energy[channel]);
        const double rmsZS = h_ohcal_channel_energy[channel]->GetRMS();
        h_CaloValid_ohcal_etaphi_ZSpedRMS->Fill(ieta, iphi, rmsZS);
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloValid::MirrorHistogram(TH1* h)
{
  const int middleBin = h->GetXaxis()->FindBin(0.0);

  for (int i = 1; i < middleBin; ++i)
  {
    const int correspondingBin = middleBin + (middleBin - i);
    const float negValue = static_cast<float>(h->GetBinContent(i));
    h->SetBinContent(correspondingBin, negValue);
  }
}

TH2* CaloValid::LogYHist2D(const std::string& name, const std::string& title, int xbins_in, double xmin, double xmax, int ybins_in, double ymin, double ymax)
{
  // note: ROOT wants Double_t*, we manage/delete manually
  const Double_t logymin = std::log10(ymin);
  const Double_t logymax = std::log10(ymax);
  const Double_t binwidth = (logymax - logymin) / static_cast<Double_t>(ybins_in);

  // allocate +2 to avoid malloc corruption in ROOT
  Double_t* ybins = new Double_t[ybins_in + 2];

  for (int i = 0; i <= ybins_in + 1; i++)
  {
    ybins[i] = std::pow(10.0, logymin + (static_cast<Double_t>(i) * binwidth));
  }

  TH2F* h = new TH2F(name.c_str(), title.c_str(),
                     xbins_in, xmin, xmax,
                     ybins_in, ybins);
  delete[] ybins;
  return h;
}

std::string CaloValid::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void CaloValid::createHistos()
{
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  h_emcal_mbd_correlation = new TH2F(std::format("{}emcal_mbd_correlation", getHistoPrefix()).c_str(), ";emcal;mbd", 100, 0, 1, 100, 0, 1);
  h_emcal_mbd_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_mbd_correlation);

  h_mbd_hits = new TH1F(std::format("{}mbd_hits", getHistoPrefix()).c_str(), "mb hits", 100, 0, 100);
  h_mbd_hits->SetDirectory(nullptr);
  hm->registerHisto(h_mbd_hits);

  h_ohcal_mbd_correlation = new TH2F(std::format("{}ohcal_mbd_correlation", getHistoPrefix()).c_str(), ";ohcal;mbd", 100, 0, 1, 100, 0, 1);
  h_ohcal_mbd_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_mbd_correlation);

  h_ihcal_mbd_correlation = new TH2F(std::format("{}ihcal_mbd_correlation", getHistoPrefix()).c_str(), ";ihcal;mbd", 100, 0, 1, 100, 0, 1);
  h_ihcal_mbd_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_mbd_correlation);

  h_emcal_hcal_correlation = new TH2F(std::format("{}emcal_hcal_correlation", getHistoPrefix()).c_str(), ";emcal;hcal", 100, 0, 1, 100, 0, 1);
  h_emcal_hcal_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_hcal_correlation);

  h_cemc_etaphi = new TH2F(std::format("{}cemc_etaphi", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256);
  h_cemc_etaphi->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi);

  h_cemc_etaphi_highhit = new TH2F(std::format("{}cemc_etaphi_highthreshold", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256);
  h_cemc_etaphi_highhit->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_highhit);

  h_ihcal_etaphi = new TH2F(std::format("{}ihcal_etaphi", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ihcal_etaphi->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi);

  h_ihcal_etaphi_highhit = new TH2F(std::format("{}ihcal_etaphi_highthreshold", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ihcal_etaphi_highhit->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_highhit);

  h_ohcal_etaphi = new TH2F(std::format("{}ohcal_etaphi", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ohcal_etaphi->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi);

  h_ohcal_etaphi_highhit = new TH2F(std::format("{}ohcal_etaphi_highthreshold", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ohcal_etaphi_highhit->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_highhit);

  h_cemc_etaphi_wQA = new TH2F(std::format("{}cemc_etaphi_wQA", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256);
  h_cemc_etaphi_wQA->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_wQA);

  h_ihcal_etaphi_wQA = new TH2F(std::format("{}ihcal_etaphi_wQA", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ihcal_etaphi_wQA->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_wQA);

  h_ohcal_etaphi_wQA = new TH2F(std::format("{}ohcal_etaphi_wQA", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ohcal_etaphi_wQA->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_wQA);

  h_ihcal_status = new TH1F(std::format("{}ihcal_status", getHistoPrefix()).c_str(), "", 256, 0, 256);
  h_ihcal_status->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_status);

  h_ohcal_status = new TH1F(std::format("{}ohcal_status", getHistoPrefix()).c_str(), "", 256, 0, 256);
  h_ohcal_status->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_status);

  h_cemc_status = new TH1F(std::format("{}cemc_status", getHistoPrefix()).c_str(), "", 256, 0, 256);
  h_cemc_status->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_status);

  h_cemc_e_chi2 = LogYHist2D(std::format("{}cemc_e_chi2", getHistoPrefix()), "", 270, -2, 25, 1000, 0.5, 4e8);
  h_cemc_e_chi2->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_e_chi2);

  h_ihcal_e_chi2 = LogYHist2D(std::format("{}ihcal_e_chi2", getHistoPrefix()), "", 270, -2, 25, 1000, 0.5, 4e8);
  h_ihcal_e_chi2->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_e_chi2);

  h_ohcal_e_chi2 = LogYHist2D(std::format("{}ohcal_e_chi2", getHistoPrefix()), "", 270, -2, 25, 1000, 0.5, 4e8);
  h_ohcal_e_chi2->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_e_chi2);

  h_cemc_etaphi_time = new TProfile2D(std::format("{}cemc_etaphi_time", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_time->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_time);

  h_cemc_etaphi_time_raw = new TProfile2D(std::format("{}cemc_etaphi_time_raw", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_time_raw->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_time_raw);

  h_cemc_etaphi_time_highhit = new TProfile2D(std::format("{}cemc_etaphi_time_highthreshold", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_time_highhit->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_time_highhit);

  h_ihcal_etaphi_time = new TProfile2D(std::format("{}ihcal_etaphi_time", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_time->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_time);

  h_ihcal_etaphi_time_raw = new TProfile2D(std::format("{}ihcal_etaphi_time_raw", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_time_raw->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_time_raw);

  h_ihcal_etaphi_time_highhit = new TProfile2D(std::format("{}ihcal_etaphi_time_highthreshold", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_time_highhit->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_time_highhit);

  h_ohcal_etaphi_time = new TProfile2D(std::format("{}ohcal_etaphi_time", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_time->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_time);

  h_ohcal_etaphi_time_raw = new TProfile2D(std::format("{}ohcal_etaphi_time_raw", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_time_raw->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_time_raw);

  h_ohcal_etaphi_time_highhit = new TProfile2D(std::format("{}ohcal_etaphi_time_highthreshold", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_time_highhit->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_time_highhit);

  h_cemc_etaphi_fracHitADC = new TProfile2D(std::format("{}cemc_etaphi_fracHitADC", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_fracHitADC->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_fracHitADC);

  h_cemc_etaphi_fracHit = new TProfile2D(std::format("{}cemc_etaphi_fracHit", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_fracHit->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_fracHit);

  h_ihcal_etaphi_fracHitADC = new TProfile2D(std::format("{}ihcal_etaphi_fracHitADC", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_fracHitADC->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_fracHitADC);

  h_ohcal_etaphi_fracHitADC = new TProfile2D(std::format("{}ohcal_etaphi_fracHitADC", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_fracHitADC->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_fracHitADC);

  h_cemc_etaphi_pedRMS = new TProfile2D(std::format("{}cemc_etaphi_pedRMS", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, 0, 1000);
  h_cemc_etaphi_pedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_pedRMS);

  h_ohcal_etaphi_pedRMS = new TProfile2D(std::format("{}ohcal_etaphi_pedRMS", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ohcal_etaphi_pedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_pedRMS);

  h_ihcal_etaphi_pedRMS = new TProfile2D(std::format("{}ihcal_etaphi_pedRMS", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ihcal_etaphi_pedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_pedRMS);

  h_cemc_etaphi_ZSpedRMS = new TProfile2D(std::format("{}cemc_etaphi_ZSpedRMS", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, 0, 1000);
  h_cemc_etaphi_ZSpedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_ZSpedRMS);

  h_ohcal_etaphi_ZSpedRMS = new TProfile2D(std::format("{}ohcal_etaphi_ZSpedRMS", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ohcal_etaphi_ZSpedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_ZSpedRMS);

  h_ihcal_etaphi_ZSpedRMS = new TProfile2D(std::format("{}ihcal_etaphi_ZSpedRMS", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ihcal_etaphi_ZSpedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_ZSpedRMS);

  h_cemc_etaphi_badChi2 = new TProfile2D(std::format("{}cemc_etaphi_badChi2", getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_badChi2->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_badChi2);

  h_ihcal_etaphi_badChi2 = new TProfile2D(std::format("{}ihcal_etaphi_badChi2", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_badChi2->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_badChi2);

  h_ohcal_etaphi_badChi2 = new TProfile2D(std::format("{}ohcal_etaphi_badChi2", getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_badChi2->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_badChi2);

  // 1D distributions
  h_InvMass = new TH1F(std::format("{}InvMass", getHistoPrefix()).c_str(), "Invariant Mass", 120, 0, 1.2);
  h_InvMass->SetDirectory(nullptr);
  hm->registerHisto(h_InvMass);

  h_channel_pedestal_0 = new TH1F(std::format("{}channel_pedestal_0", getHistoPrefix()).c_str(), "Test Pedestal", 1000, -500., 500.);
  h_channel_pedestal_0->SetDirectory(nullptr);
  hm->registerHisto(h_channel_pedestal_0);

  // vertex distributions
  h_vtx_z_raw = new TH1D(std::format("{}vtx_z_raw", getHistoPrefix()).c_str(), "hvtx_z_raw", 201, -100.5, 100.5);
  h_vtx_z_raw->SetDirectory(nullptr);
  hm->registerHisto(h_vtx_z_raw);

  h_vtx_z_cut = new TH1D(std::format("{}vtx_z_cut", getHistoPrefix()).c_str(), "hvtx_z_cut", 201, -100.5, 100.5);
  h_vtx_z_cut->SetDirectory(nullptr);
  hm->registerHisto(h_vtx_z_cut);

  // raw timing information
  h_emcaltime_cut = new TH1D(std::format("{}emcaltime_cut", getHistoPrefix()).c_str(), "hemcaltime_cut", 50, -17.5, 32.5);
  h_emcaltime_cut->SetDirectory(nullptr);
  hm->registerHisto(h_emcaltime_cut);

  h_ihcaltime_cut = new TH1D(std::format("{}ihcaltime_cut", getHistoPrefix()).c_str(), "hihcaltime_cut", 50, -17.5, 32.5);
  h_ihcaltime_cut->SetDirectory(nullptr);
  hm->registerHisto(h_ihcaltime_cut);

  h_ohcaltime_cut = new TH1D(std::format("{}ohcaltime_cut", getHistoPrefix()).c_str(), "hohcaltime_cut", 50, -17.5, 32.5);
  h_ohcaltime_cut->SetDirectory(nullptr);
  hm->registerHisto(h_ohcaltime_cut);

  // extracted timing information
  h_emcaltime = new TH1D(std::format("{}emcaltime", getHistoPrefix()).c_str(), "hemcaltime", 50, -17.5, 32.5);
  h_emcaltime->SetDirectory(nullptr);
  hm->registerHisto(h_emcaltime);

  h_ihcaltime = new TH1D(std::format("{}ihcaltime", getHistoPrefix()).c_str(), "hihcaltime", 50, -17.5, 32.5);
  h_ihcaltime->SetDirectory(nullptr);
  hm->registerHisto(h_ihcaltime);

  h_ohcaltime = new TH1D(std::format("{}ohcaltime", getHistoPrefix()).c_str(), "hohcaltime", 50, -17.5, 32.5);
  h_ohcaltime->SetDirectory(nullptr);
  hm->registerHisto(h_ohcaltime);

  h_emcal_tower_e = new TH1F(std::format("{}emcal_tower_e", getHistoPrefix()).c_str(), "emcal_tower_e", 5000, -0.1, 1);
  h_emcal_tower_e->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_tower_e);

  h_emcal_tower_e_wide_range = new TH1F(std::format("{}emcal_tower_e_wide_range", getHistoPrefix()).c_str(), "emcal_tower_e_wide_range", 1000, -10., 40.);
  h_emcal_tower_e_wide_range->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_tower_e_wide_range);

  h_emcal_tower_e_saturated = new TH1F(std::format("{}emcal_tower_e_saturated", getHistoPrefix()).c_str(), "emcal_tower_e_saturated", 1000, -10., 40.);
  h_emcal_tower_e_saturated->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_tower_e_saturated);

  h_ihcal_tower_e = new TH1F(std::format("{}ihcal_tower_e", getHistoPrefix()).c_str(), "ihcal_tower_e", 5000, -0.1, 1);
  h_ihcal_tower_e->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_tower_e);

  h_ihcal_tower_e_wide_range = new TH1F(std::format("{}ihcal_tower_e_wide_range", getHistoPrefix()).c_str(), "ihcal_tower_e_wide_range", 1000, -10., 40.);
  h_ihcal_tower_e_wide_range->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_tower_e_wide_range);

  h_ihcal_tower_e_saturated = new TH1F(std::format("{}ihcal_tower_e_saturated", getHistoPrefix()).c_str(), "ihcal_tower_e_saturated", 1000, -10., 40.);
  h_ihcal_tower_e_saturated->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_tower_e_saturated);

  h_ohcal_tower_e = new TH1F(std::format("{}ohcal_tower_e", getHistoPrefix()).c_str(), "ohcal_tower_e", 5000, -0.1, 1);
  h_ohcal_tower_e->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_tower_e);

  h_ohcal_tower_e_wide_range = new TH1F(std::format("{}ohcal_tower_e_wide_range", getHistoPrefix()).c_str(), "ohcal_tower_e_wide_range", 1000, -10., 40.);
  h_ohcal_tower_e_wide_range->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_tower_e_wide_range);

  h_ohcal_tower_e_saturated = new TH1F(std::format("{}ohcal_tower_e_saturated", getHistoPrefix()).c_str(), "ohcal_tower_e_saturated", 1000, -10., 40.);
  h_ohcal_tower_e_saturated->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_tower_e_saturated);

  h_etaphi_clus = new TH2F(std::format("{}etaphi_clus", getHistoPrefix()).c_str(), "", 140, -1.2, 1.2, 64, -1 * M_PI, M_PI);
  h_etaphi_clus->SetDirectory(nullptr);
  hm->registerHisto(h_etaphi_clus);

  h_clusE = new TH1F(std::format("{}clusE", getHistoPrefix()).c_str(), "", 100, 0, 10);
  h_clusE->SetDirectory(nullptr);
  hm->registerHisto(h_clusE);

  h_triggerVec = new TH1F(std::format("{}triggerVec", getHistoPrefix()).c_str(), "", 64, 0, 64);
  h_triggerVec->SetDirectory(nullptr);
  hm->registerHisto(h_triggerVec);

  //---------EMCal channel histos--------//
  {
    const int size = 128 * 192;
    for (int channel = 0; channel < size; channel++)
    {
      std::string hname = std::format("h_cemc_channel_pedestal_{}", channel);
      h_cemc_channel_pedestal[channel] = new TH1F(hname.c_str(), hname.c_str(), 2000, -0.5, 2000.5);
      h_cemc_channel_pedestal[channel]->SetDirectory(nullptr);

      std::string hnameE = std::format("h_cemc_channel_energy_{}", channel);
      h_cemc_channel_energy[channel] = new TH1F(hnameE.c_str(), hnameE.c_str(), 1000, -50, 50);
      h_cemc_channel_energy[channel]->SetDirectory(nullptr);
    }
  }
  //--------OHCal channel histos--------//
  {
    const int size = 32 * 48;
    for (int channel = 0; channel < size; channel++)
    {
      std::string hname = std::format("h_ohcal_channel_pedestal_{}", channel);
      h_ohcal_channel_pedestal[channel] = new TH1F(hname.c_str(), hname.c_str(), 2000, -0.5, 2000.5);
      h_ohcal_channel_pedestal[channel]->SetDirectory(nullptr);

      std::string hnameE = std::format("h_ohcal_channel_energy_{}", channel);
      h_ohcal_channel_energy[channel] = new TH1F(hnameE.c_str(), hnameE.c_str(), 1000, -50, 50);
      h_ohcal_channel_energy[channel]->SetDirectory(nullptr);
    }
  }
  //--------IHCal channel histos-------//
  {
    const int size = 32 * 48;
    for (int channel = 0; channel < size; channel++)
    {
      std::string hname = std::format("h_ihcal_channel_pedestal_{}", channel);
      h_ihcal_channel_pedestal[channel] = new TH1F(hname.c_str(), hname.c_str(), 2000, -0.5, 2000.5);
      h_ihcal_channel_pedestal[channel]->SetDirectory(nullptr);

      std::string hnameE = std::format("h_ihcal_channel_energy_{}", channel);
      h_ihcal_channel_energy[channel] = new TH1F(hnameE.c_str(), hnameE.c_str(), 1000, -50, 50);
      h_ihcal_channel_energy[channel]->SetDirectory(nullptr);
    }
  }

  h_pi0_trigIB_mass = new TH3F(
      "h_pi0_trigIB_mass",
      ";Trigger Bit; iIB (eta*32 + phi); #pi^{0} Mass (GeV/c^{2})",
      64, -0.5, 63.5,    // triggers
      384, -0.5, 383.5,  // IB index
      120, 0.0, 1.2      // mass range
  );
  h_pi0_trigIB_mass->SetDirectory(nullptr);
  hm->registerHisto(h_pi0_trigIB_mass);

  // Trigger QA
  pr_ldClus_trig = new TProfile("pr_CaloValid_ldClus_trig", "", 64, 0, 64, 0, 10);
  for (int i = 0; i < 64; i++)
  {
    if (!inTrigOfInterest(trigOfInterest, i))
    {
      continue;
    }

    h_edist[i] = new TH2F(
        std::format("h_CaloValid_edist_trig{}", i).c_str(), "",
        64, -1.2, 1.2, 128, -3.1415, 3.1415);

    h_ldClus_trig[i] = new TH1F(
        std::format("h_CaloValid_ldClus_trig{}", i).c_str(), "",
        18, 1, 10);

    pr_evtNum_ldClus_trig[i] = new TProfile(
        std::format("pr_CaloValid_evtNum_ldClus_trig{}", i).c_str(),
        "", 100000, 0, 100000, 0, 10);

    pr_rejection[i] = new TProfile(
        std::format("pr_CaloValid_rejection_trig{}", i).c_str(),
        "", 100000, 0, 100000, 0, 50000);

    pr_livetime[i] = new TProfile(
        std::format("pr_CaloValid_livetime_trig{}", i).c_str(),
        "", 100000, 0, 100000, 0, 10);

    hm->registerHisto(h_edist[i]);
    hm->registerHisto(h_ldClus_trig[i]);
    hm->registerHisto(pr_evtNum_ldClus_trig[i]);
    hm->registerHisto(pr_rejection[i]);
    hm->registerHisto(pr_livetime[i]);
  }

  hm->registerHisto(h_triggerVec);
  hm->registerHisto(pr_ldClus_trig);
}
