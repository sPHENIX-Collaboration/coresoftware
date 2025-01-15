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

#include <ffarawobjects/Gl1Packet.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE


#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <cassert>
#include <cmath>  // for log10, pow, sqrt, abs, M_PI
#include <cstdint>
#include <iostream>  // for operator<<, endl, basic_...
#include <limits>
#include <map>  // for operator!=, _Rb_tree_con...
#include <string>
#include <utility>  // for pair

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

int CaloValid::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;
  //  std::cout << "In process_event" << std::endl;
  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::process_towers(PHCompositeNode* topNode)
{
  //---------------------------Event header--------------------------------//
  EventHeader* eventheader =
      findNode::getClass<EventHeader>(topNode, "EventHeader");
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
  //  std::cout << "In process_towers" << std::endl;
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  float totalcemc = 0.;
  float totalihcal = 0.;
  float totalohcal = 0.;
  float totalmbd = 0.;

  float emcaldownscale = 100000. / 800.;
  float ihcaldownscale = 4000. / 300.;
  float ohcaldownscale = 25000. / 600.;
  float mbddownscale = 200.0;

  float adc_threshold = 100.;

  float emcal_hit_threshold = 0.5;  // GeV
  float ohcal_hit_threshold = 0.5;
  float ihcal_hit_threshold = 0.25;

  //----------------------------------vertex------------------------------------------------------//
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
      vtx_z = vtx->get_z();
    }
    h_vtx_z_raw->Fill(vtx_z);
  }

  //--------------------------- trigger and GL1-------------------------------//
  if (trigAna)
  {
    trigAna->decodeTriggers(topNode);
  }
  else
  {
    if (m_debug)
    {
    std::cout << "[ERROR] No TriggerAnalyzer pointer!\n";
    }
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
  long long int raw[64] = {0};
  long long int live[64] = {0};
  // long long int scaled[64] = { 0 };
  Gl1Packet* gl1PacketInfo =
      findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1PacketInfo)
  {
    std::cout << PHWHERE << "GlobalQA::process_event: GL1Packet node is missing"
              << std::endl;
  }
  uint64_t triggervec = 0;
  if (gl1PacketInfo)
  {
    triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++)
    {
      bool trig_decision = ((triggervec & 0x1U) == 0x1U);
      scaledBits[i] = trig_decision;

      raw[i] = gl1PacketInfo->lValue(i, 0);
      live[i] = gl1PacketInfo->lValue(i, 1);
      // scaled[i] = gl1PacketInfo->lValue(i, 2);

      if (trig_decision)
      {
        h_triggerVec->Fill(i);
      }
      triggervec = (triggervec >> 1U) & 0xffffffffU;
    }
    //triggervec = gl1PacketInfo->getScaledVector(); commented out to get rid of never used warning using clang -tidy
  }

  //---------------------------calibrated towers-------------------------------//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        float offlineenergy = tower->get_energy();
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        int _time = tower->get_time();
        h_cemc_e_chi2->Fill(offlineenergy, tower->get_chi2());
        float _timef = tower->get_time_float();
        h_emcaltime_cut->Fill(_time);
        bool isGood = tower->get_isGood();
        uint8_t status = tower->get_status();
        h_emcal_tower_e->Fill(offlineenergy);
        float pedestal = tower->get_pedestal();
        h_cemc_channel_pedestal[channel]->Fill(pedestal);

        for (int is = 0; is < 8; is++)
        {
          if (status & 1U)  // clang-tidy mark 1 as unsigned
          {
            h_cemc_status->Fill(is);
          }
          status = status >> 1U;  // clang-tidy mark 1 as unsigned
        }

        totalcemc += offlineenergy;
        h_emcaltime->Fill(_time);
        if (offlineenergy > emcal_hit_threshold)
        {
          h_cemc_etaphi_time->Fill(ieta, iphi, _timef);
          h_cemc_etaphi->Fill(ieta, iphi);
          if (isGood && (scaledBits[10] || scaledBits[11]))
          {
            h_cemc_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
          }
          if (tower->get_isBadChi2())
          {
            h_cemc_etaphi_badChi2->Fill(ieta, iphi, 1);
          }
          else
          {
            h_cemc_etaphi_badChi2->Fill(ieta, iphi, 0);
          }
        }
        if (offlineenergy > 0.25)
        {
          h_cemc_etaphi_fracHit->Fill(ieta, iphi, 1);
        }
        else
        {
          h_cemc_etaphi_fracHit->Fill(ieta, iphi, 0);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        float offlineenergy = tower->get_energy();
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        int _time = tower->get_time();
        float _timef = tower->get_time_float();
        h_ihcaltime_cut->Fill(_time);
        h_ihcal_e_chi2->Fill(offlineenergy, tower->get_chi2());
        bool isGood = tower->get_isGood();
        h_ihcal_status->Fill(tower->get_status());
        float pedestal = tower->get_pedestal();
        h_ihcal_channel_pedestal[channel]->Fill(pedestal);

        uint8_t status = tower->get_status();
        for (int is = 0; is < 8; is++)
        {
          if (status & 1U)  // clang-tidy mark 1 as unsigned
          {
            h_ihcal_status->Fill(is);
          }
          status = status >> 1U;  // clang-tidy mark 1 as unsigned
        }

        totalihcal += offlineenergy;
        h_ihcaltime->Fill(_time);

        if (offlineenergy > ihcal_hit_threshold)
        {
          h_ihcal_etaphi->Fill(ieta, iphi);
          h_ihcal_etaphi_time->Fill(ieta, iphi, _timef);
          if (isGood && (scaledBits[10] || scaledBits[11]))
          {
            h_ihcal_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
          }
          if (tower->get_isBadChi2())
          {
            h_ihcal_etaphi_badChi2->Fill(ieta, iphi, 1);
          }
          else
          {
            h_ihcal_etaphi_badChi2->Fill(ieta, iphi, 0);
          }
        }
      }
    }
  }

  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        float offlineenergy = tower->get_energy();
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        int _time = tower->get_time();
        float _timef = tower->get_time_float();
        h_ohcaltime_cut->Fill(_time);
        h_ohcal_e_chi2->Fill(offlineenergy, tower->get_chi2());
        bool isGood = tower->get_isGood();
        h_ohcal_status->Fill(tower->get_status());
        float pedestal = tower->get_pedestal();
        h_ohcal_channel_pedestal[channel]->Fill(pedestal);

        uint8_t status = tower->get_status();
        for (int is = 0; is < 8; is++)
        {
          if (status & 1U)  // clang-tidy mark 1 as unsigned
          {
            h_ohcal_status->Fill(is);
          }
          status = status >> 1U;  // clang-tidy mark 1 as unsigned
        }

        totalohcal += offlineenergy;
        h_ohcaltime->Fill(_time);

        if (offlineenergy > ohcal_hit_threshold)
        {
          h_ohcal_etaphi->Fill(ieta, iphi);
          h_ohcal_etaphi_time->Fill(ieta, iphi, _timef);
          if (isGood && (scaledBits[10] || scaledBits[11]))
          {
            h_ohcal_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
          }
          if (tower->get_isBadChi2())
          {
            h_ohcal_etaphi_badChi2->Fill(ieta, iphi, 1);
          }
          else
          {
            h_ohcal_etaphi_badChi2->Fill(ieta, iphi, 0);
          }
        }
      }
    }
  }

  //-------------------------- raw tower ------------------------------//

  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float raw_time = tower->get_time_float();
        if (tower->get_isZS())
        {
          h_cemc_channel_energy[channel]->Fill(tower->get_energy());
        }

        float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_cemc_etaphi_fracHitADC->Fill(ieta, iphi, 1);
          h_cemc_etaphi_time_raw->Fill(ieta, iphi, raw_time);
        }
        else
        {
          h_cemc_etaphi_fracHitADC->Fill(ieta, iphi, 0);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALOUT");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float raw_time = tower->get_time_float();
        if (tower->get_isZS())
        {
          h_ohcal_channel_energy[channel]->Fill(tower->get_energy());
        }

        float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_ohcal_etaphi_time_raw->Fill(ieta, iphi, raw_time);
          h_ohcal_etaphi_fracHitADC->Fill(ieta, iphi, 1);
        }
        else
        {
          h_ohcal_etaphi_fracHitADC->Fill(ieta, iphi, 0);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        float raw_time = tower->get_time_float();
        int iphi = towers->getTowerPhiBin(towerkey);
        if (tower->get_isZS())
        {
          h_ihcal_channel_energy[channel]->Fill(tower->get_energy());
        }

        float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_ihcal_etaphi_time_raw->Fill(ieta, iphi, raw_time);
          h_ihcal_etaphi_fracHitADC->Fill(ieta, iphi, 1);
        }
        else
        {
          h_ihcal_etaphi_fracHitADC->Fill(ieta, iphi, 0);
        }
      }
    }
  }

  //--------------------------- MBD ----------------------------------------//
  MbdPmtContainer* bbcpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!bbcpmts)
  {
    std::cout << "CaloVald::process_event: Could not find MbdPmtContainer," << std::endl;
    // return Fun4AllReturnCodes::ABORTEVENT;
  }

  int hits = 0;
  if (bbcpmts)
  {
    int nPMTs = bbcpmts->get_npmt();
    for (int i = 0; i < nPMTs; i++)
    {
      MbdPmtHit* mbdpmt = bbcpmts->get_pmt(i);
      float pmtadc = mbdpmt->get_q();
      totalmbd += pmtadc;
      if (pmtadc > 0.4)
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

  //------------------------------ clusters & pi0 ------------------------------//
  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "CaloValid::funkyCaloStuff::process_event - Fatal Error - CLUSTER_CEMC node is missing. " << std::endl;
    return 0;
  }

  //////////////////////////////////////////
  // geometry for hot tower/cluster masking
  std::string towergeomnodename = "TOWERGEOM_CEMC";
  RawTowerGeomContainer* m_geometry = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!m_geometry)
  {
    std::cout << Name() << "::"
              << "CreateNodeTree"
              << ": Could not find node " << towergeomnodename << std::endl;
    gSystem->Exit(1);
  }

  // cuts
  float emcMinClusE1 = 1.3;  // 0.5;
  float emcMinClusE2 = 0.7;  // 0.5;
  float emcMaxClusE = 100;
  float maxAlpha = 0.6;

  if (totalcemc < 0.2 * emcaldownscale)
  {
    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter;
    RawClusterContainer::ConstIterator clusterIter2;

    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
    {
      RawCluster* recoCluster = clusterIter->second;

      CLHEP::Hep3Vector vertex(0, 0, 0);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

      float clusE = E_vec_cluster.mag();
      float clus_eta = E_vec_cluster.pseudoRapidity();
      float clus_phi = E_vec_cluster.phi();
      float clus_pt = E_vec_cluster.perp();
      float clus_chisq = recoCluster->get_chi2();

      h_clusE->Fill(clusE);

      if (clusE < emcMinClusE1 || clusE > emcMaxClusE)
      {
        continue;
      }
      if (clus_chisq > 4)
      {
        continue;
      }

      h_etaphi_clus->Fill(clus_eta, clus_phi);

      TLorentzVector photon1;
      photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);

      for (clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; clusterIter2++)
      {
        if (clusterIter == clusterIter2)
        {
          continue;
        }
        RawCluster* recoCluster2 = clusterIter2->second;

        CLHEP::Hep3Vector vertex2(0, 0, 0);
        CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*recoCluster2, vertex2);

        float clus2E = E_vec_cluster2.mag();
        float clus2_eta = E_vec_cluster2.pseudoRapidity();
        float clus2_phi = E_vec_cluster2.phi();
        float clus2_pt = E_vec_cluster2.perp();
        float clus2_chisq = recoCluster2->get_chi2();

        if (clus2E < emcMinClusE2 || clus2E > emcMaxClusE)
        {
          continue;
        }
        if (clus2_chisq > 4)
        {
          continue;
        }

        TLorentzVector photon2;
        photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);

        if (sqrt(pow(clusE - clus2E, 2)) / (clusE + clus2E) > maxAlpha)
        {
          continue;
        }

        TLorentzVector pi0 = photon1 + photon2;
        float pi0Mass = pi0.M();
        unsigned int lt_eta = recoCluster->get_lead_tower().first;
        unsigned int lt_phi = recoCluster->get_lead_tower().second;

        int ld_ib_eta = lt_eta / 8;
        int ld_ib_phi = lt_phi / 8;
        int IB_num    = ld_ib_eta * 12 + ld_ib_phi;

        for (int bit : scaledActiveBits)
        {
          if (std::find(triggerIndices.begin(), triggerIndices.end(), bit)
              == triggerIndices.end())
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
    } // end cluster loop
  }

  //----------------- Trigger / alignment ----------------------------//
  float leading_cluster_ecore = 0;
  float leading_cluster_eta = 0;
  float leading_cluster_phi = 0;
  int evtNum_overK = event_number / 1000;

  if (clusterContainer)
  {
    RawClusterContainer::ConstRange clusterEnd =
        clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter;
    RawClusterContainer::ConstIterator clusterIter2;

    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second;
         clusterIter++)
    {
      RawCluster* recoCluster = clusterIter->second;
      if (recoCluster->get_chi2() > 2)
      {
        continue;
      }

      CLHEP::Hep3Vector vertex(0, 0, 0);
      CLHEP::Hep3Vector E_vec_cluster =
          RawClusterUtility::GetECoreVec(*recoCluster, vertex);

      float clusE = E_vec_cluster.mag();
      float clusEta = E_vec_cluster.pseudoRapidity();
      float clusPhi = E_vec_cluster.phi();
      if (clusE > leading_cluster_ecore)
      {
        leading_cluster_ecore = clusE;
        leading_cluster_eta = clusEta;
        leading_cluster_phi = clusPhi;
      }
    }
    for (int i = 0; i < 64; i++)
    {
      if (scaledBits[i])
      {
        pr_ldClus_trig->Fill(i, leading_cluster_ecore);
        if (!(std::find(trigOfInterest.begin(), trigOfInterest.end(), i) != trigOfInterest.end()))
        {
          continue;
        }
        h_edist[i]->Fill(leading_cluster_eta, leading_cluster_phi);
        h_ldClus_trig[i]->Fill(leading_cluster_ecore);
        pr_evtNum_ldClus_trig[i]->Fill(evtNum_overK, leading_cluster_ecore);
        if (raw[i] > 0)
        {
          pr_rejection[i]->Fill(evtNum_overK,
                                (float) raw[10] / (float) raw[i]);
          pr_livetime[i]->Fill(evtNum_overK,
                               (float) live[i] / (float) raw[i]);
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::End(PHCompositeNode* topNode)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  //------EmCal-----//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if (towers)
    {
      int size = towers->size();

      auto h_CaloValid_cemc_etaphi_pedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%scemc_etaphi_pedRMS") % getHistoPrefix()).c_str()));
      auto h_CaloValid_cemc_etaphi_ZSpedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%scemc_etaphi_ZSpedRMS") % getHistoPrefix()).c_str()));

      for (int channel = 0; channel < size; channel++)
      {
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float ped_rms = h_cemc_channel_pedestal[channel]->GetRMS();
        h_CaloValid_cemc_etaphi_pedRMS->Fill(ieta, iphi, ped_rms);
        MirrorHistogram(h_cemc_channel_energy[channel]);
        double rmsZS = h_cemc_channel_energy[channel]->GetRMS();
        h_CaloValid_cemc_etaphi_ZSpedRMS->Fill(ieta, iphi, rmsZS);
      }
    }
  }
  //------IHCal------//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if (towers)
    {
      int size = towers->size();

      auto h_CaloValid_ihcal_etaphi_pedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sihcal_etaphi_pedRMS") % getHistoPrefix()).c_str()));
      auto h_CaloValid_ihcal_etaphi_ZSpedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sihcal_etaphi_ZSpedRMS") % getHistoPrefix()).c_str()));

      for (int channel = 0; channel < size; channel++)
      {
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float ped_rms = h_ihcal_channel_pedestal[channel]->GetRMS();
        h_CaloValid_ihcal_etaphi_pedRMS->Fill(ieta, iphi, ped_rms);
        MirrorHistogram(h_ihcal_channel_energy[channel]);
        double rmsZS = h_ihcal_channel_energy[channel]->GetRMS();
        h_CaloValid_ihcal_etaphi_ZSpedRMS->Fill(ieta, iphi, rmsZS);
      }
    }
  }
  //------OHCal-----//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if (towers)
    {
      int size = towers->size();

      auto h_CaloValid_ohcal_etaphi_pedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sohcal_etaphi_pedRMS") % getHistoPrefix()).c_str()));
      auto h_CaloValid_ohcal_etaphi_ZSpedRMS = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%oohcal_etaphi_ZSpedRMS") % getHistoPrefix()).c_str()));

      for (int channel = 0; channel < size; channel++)
      {
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float ped_rms = h_ohcal_channel_pedestal[channel]->GetRMS();
        h_CaloValid_ohcal_etaphi_pedRMS->Fill(ieta, iphi, ped_rms);
        MirrorHistogram(h_ohcal_channel_energy[channel]);
        double rmsZS = h_ohcal_channel_energy[channel]->GetRMS();
        h_CaloValid_ohcal_etaphi_ZSpedRMS->Fill(ieta, iphi, rmsZS);
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloValid::MirrorHistogram(TH1* h)
{
  int middleBin = h->GetXaxis()->FindBin(0.0);

  for (int i = 1; i < middleBin; ++i)
  {
    int correspondingBin = middleBin + (middleBin - i);
    float negValue = h->GetBinContent(i);
    h->SetBinContent(correspondingBin, negValue);
  }
}

TH2* CaloValid::LogYHist2D(const std::string& name, const std::string& title, int xbins_in, double xmin, double xmax, int ybins_in, double ymin, double ymax)
{
  Double_t logymin = std::log10(ymin);
  Double_t logymax = std::log10(ymax);
  Double_t binwidth = (logymax - logymin) / ybins_in;
  Double_t *ybins = new Double_t[ybins_in + 2]; //allocate 1 extra bin to fix "malloc()L memory corruption crash

  for (Int_t i = 0; i <= ybins_in + 1; i++)
  {
    ybins[i] = pow(10, logymin + i * binwidth);
  }

  TH2F* h = new TH2F(name.c_str(), title.c_str(), xbins_in, xmin, xmax, ybins_in, ybins);
  delete [] ybins;
  return h;
}
std::string CaloValid::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void CaloValid::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create and register your histos (all types) here
  h_emcal_mbd_correlation = new TH2F(boost::str(boost::format("%semcal_mbd_correlation") % getHistoPrefix()).c_str(), ";emcal;mbd", 100, 0, 1, 100, 0, 1);
  h_emcal_mbd_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_mbd_correlation);

  h_mbd_hits = new TH1F(boost::str(boost::format("%smbd_hits") % getHistoPrefix()).c_str(), "mb hits", 100, 0, 100);
  h_mbd_hits->SetDirectory(nullptr);
  hm->registerHisto(h_mbd_hits);

  h_ohcal_mbd_correlation = new TH2F(boost::str(boost::format("%sohcal_mbd_correlation") % getHistoPrefix()).c_str(), ";ohcal;mbd", 100, 0, 1, 100, 0, 1);
  h_ohcal_mbd_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_mbd_correlation);

  h_ihcal_mbd_correlation = new TH2F(boost::str(boost::format("%sihcal_mbd_correlation") % getHistoPrefix()).c_str(), ";ihcal;mbd", 100, 0, 1, 100, 0, 1);
  h_ihcal_mbd_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_mbd_correlation);

  h_emcal_hcal_correlation = new TH2F(boost::str(boost::format("%semcal_hcal_correlation") % getHistoPrefix()).c_str(), ";emcal;hcal", 100, 0, 1, 100, 0, 1);
  h_emcal_hcal_correlation->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_hcal_correlation);

  h_cemc_etaphi = new TH2F(boost::str(boost::format("%scemc_etaphi") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256);
  h_cemc_etaphi->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi);

  h_ihcal_etaphi = new TH2F(boost::str(boost::format("%sihcal_etaphi") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ihcal_etaphi->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi);

  h_ohcal_etaphi = new TH2F(boost::str(boost::format("%sohcal_etaphi") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ohcal_etaphi->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi);

  h_cemc_etaphi_wQA = new TH2F(boost::str(boost::format("%scemc_etaphi_wQA") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256);
  h_cemc_etaphi_wQA->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_wQA);

  h_ihcal_etaphi_wQA = new TH2F(boost::str(boost::format("%sihcal_etaphi_wQA") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ihcal_etaphi_wQA->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_wQA);

  h_ohcal_etaphi_wQA = new TH2F(boost::str(boost::format("%sohcal_etaphi_wQA") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
  h_ohcal_etaphi_wQA->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_wQA);

  h_ihcal_status = new TH1F(boost::str(boost::format("%sihcal_status") % getHistoPrefix()).c_str(), "", 256, 0, 256);
  h_ihcal_status->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_status);

  h_ohcal_status = new TH1F(boost::str(boost::format("%sohcal_status") % getHistoPrefix()).c_str(), "", 256, 0, 256);
  h_ohcal_status->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_status);

  h_cemc_status = new TH1F(boost::str(boost::format("%scemc_status") % getHistoPrefix()).c_str(), "", 256, 0, 256);
  h_cemc_status->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_status);

  h_cemc_e_chi2 = LogYHist2D(boost::str(boost::format("%scemc_e_chi2") % getHistoPrefix()).c_str(), "", 270, -2, 25, 1000, 0.5, 4e8);
  h_cemc_e_chi2->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_e_chi2);

  h_ihcal_e_chi2 = LogYHist2D(boost::str(boost::format("%sihcal_e_chi2") % getHistoPrefix()).c_str(), "", 270, -2, 25, 1000, 0.5, 4e8);
  h_ihcal_e_chi2->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_e_chi2);

  h_ohcal_e_chi2 = LogYHist2D(boost::str(boost::format("%sohcal_e_chi2") % getHistoPrefix()).c_str(), "", 270, -2, 25, 1000, 0.5, 4e8);
  h_ohcal_e_chi2->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_e_chi2);

  h_cemc_etaphi_time = new TProfile2D(boost::str(boost::format("%scemc_etaphi_time") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_time->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_time);

  h_cemc_etaphi_time_raw = new TProfile2D(boost::str(boost::format("%scemc_etaphi_time_raw") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_time_raw->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_time_raw);

  h_ihcal_etaphi_time = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_time") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_time->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_time);

  h_ihcal_etaphi_time_raw = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_time_raw") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_time_raw->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_time_raw);

  h_ohcal_etaphi_time = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_time") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_time->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_time);

  h_ohcal_etaphi_time_raw = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_time_raw") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_time_raw->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_time_raw);

  h_cemc_etaphi_fracHitADC = new TProfile2D(boost::str(boost::format("%scemc_etaphi_fracHitADC") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_fracHitADC->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_fracHitADC);

  h_cemc_etaphi_fracHit = new TProfile2D(boost::str(boost::format("%scemc_etaphi_fracHit") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_fracHit->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_fracHit);

  h_ihcal_etaphi_fracHitADC = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_fracHitADC") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_fracHitADC->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_fracHitADC);

  h_ohcal_etaphi_fracHitADC = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_fracHitADC") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_fracHitADC->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_fracHitADC);

  h_cemc_etaphi_pedRMS = new TProfile2D(boost::str(boost::format("%scemc_etaphi_pedRMS") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, 0, 1000);
  h_cemc_etaphi_pedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_pedRMS);

  h_ohcal_etaphi_pedRMS = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_pedRMS") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ohcal_etaphi_pedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_pedRMS);

  h_ihcal_etaphi_pedRMS = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_pedRMS") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ihcal_etaphi_pedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_pedRMS);

  h_cemc_etaphi_ZSpedRMS = new TProfile2D(boost::str(boost::format("%scemc_etaphi_ZSpedRMS") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, 0, 1000);
  h_cemc_etaphi_ZSpedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_ZSpedRMS);

  h_ohcal_etaphi_ZSpedRMS = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_ZSpedRMS") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ohcal_etaphi_ZSpedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_ZSpedRMS);

  h_ihcal_etaphi_ZSpedRMS = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_ZSpedRMS") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 1000);
  h_ihcal_etaphi_ZSpedRMS->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_ZSpedRMS);

  h_cemc_etaphi_badChi2 = new TProfile2D(boost::str(boost::format("%scemc_etaphi_badChi2") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_badChi2->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_badChi2);

  h_ihcal_etaphi_badChi2 = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_badChi2") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_badChi2->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_badChi2);

  h_ohcal_etaphi_badChi2 = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_badChi2") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_badChi2->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_badChi2);

  // 1D distributions
  h_InvMass = new TH1F(boost::str(boost::format("%sInvMass") % getHistoPrefix()).c_str(), "Invariant Mass", 120, 0, 1.2);
  h_InvMass->SetDirectory(nullptr);
  hm->registerHisto(h_InvMass);

  // for (int channel = 0; channel < 128*192; channel++) {
  h_channel_pedestal_0 = new TH1F(boost::str(boost::format("%schannel_pedestal_0") % getHistoPrefix()).c_str(), "Test Pedestal", 1000, -500., 500.);
  h_channel_pedestal_0->SetDirectory(nullptr);
  hm->registerHisto(h_channel_pedestal_0);

  // vertex distributions
  h_vtx_z_raw = new TH1D(boost::str(boost::format("%svtx_z_raw") % getHistoPrefix()).c_str(), "hvtx_z_raw", 201, -100.5, 100.5);
  h_vtx_z_raw->SetDirectory(nullptr);
  hm->registerHisto(h_vtx_z_raw);

  h_vtx_z_cut = new TH1D(boost::str(boost::format("%svtx_z_cut") % getHistoPrefix()).c_str(), "hvtx_z_cut", 201, -100.5, 100.5);
  h_vtx_z_cut->SetDirectory(nullptr);
  hm->registerHisto(h_vtx_z_cut);

  // raw timing information
  h_emcaltime_cut = new TH1D(boost::str(boost::format("%semcaltime_cut") % getHistoPrefix()).c_str(), "hemcaltime_cut", 50, -17.5, 32.5);
  h_emcaltime_cut->SetDirectory(nullptr);
  hm->registerHisto(h_emcaltime_cut);

  h_ihcaltime_cut = new TH1D(boost::str(boost::format("%sihcaltime_cut") % getHistoPrefix()).c_str(), "hihcaltime_cut", 50, -17.5, 32.5);
  h_ihcaltime_cut->SetDirectory(nullptr);
  hm->registerHisto(h_ihcaltime_cut);

  h_ohcaltime_cut = new TH1D(boost::str(boost::format("%sohcaltime_cut") % getHistoPrefix()).c_str(), "hohcaltime_cut", 50, -17.5, 32.5);
  h_ohcaltime_cut->SetDirectory(nullptr);
  hm->registerHisto(h_ohcaltime_cut);

  // extracted timing information
  h_emcaltime = new TH1D(boost::str(boost::format("%semcaltime") % getHistoPrefix()).c_str(), "hemcaltime", 50, -17.5, 32.5);
  h_emcaltime->SetDirectory(nullptr);
  hm->registerHisto(h_emcaltime);

  h_ihcaltime = new TH1D(boost::str(boost::format("%sihcaltime") % getHistoPrefix()).c_str(), "hihcaltime", 50, -17.5, 32.5);
  h_ihcaltime->SetDirectory(nullptr);
  hm->registerHisto(h_ihcaltime);

  h_ohcaltime = new TH1D(boost::str(boost::format("%sohcaltime") % getHistoPrefix()).c_str(), "hohcaltime", 50, -17.5, 32.5);
  h_ohcaltime->SetDirectory(nullptr);
  hm->registerHisto(h_ohcaltime);

  h_emcal_tower_e = new TH1F(boost::str(boost::format("%semcal_tower_e") % getHistoPrefix()).c_str(), "emcal_tower_e", 5000, -0.1, 1);
  h_emcal_tower_e->SetDirectory(nullptr);
  hm->registerHisto(h_emcal_tower_e);

  // cluster QA
  h_etaphi_clus = new TH2F(boost::str(boost::format("%setaphi_clus") % getHistoPrefix()).c_str(), "", 140, -1.2, 1.2, 64, -1 * M_PI, M_PI);
  h_etaphi_clus->SetDirectory(nullptr);
  hm->registerHisto(h_etaphi_clus);

  h_clusE = new TH1F(boost::str(boost::format("%sclusE") % getHistoPrefix()).c_str(), "", 100, 0, 10);
  h_clusE->SetDirectory(nullptr);
  hm->registerHisto(h_clusE);

  h_triggerVec = new TH1F(boost::str(boost::format("%striggerVec") % getHistoPrefix()).c_str(), "", 64, 0, 64);
  h_triggerVec->SetDirectory(nullptr);
  hm->registerHisto(h_triggerVec);

  //---------EMCal--------//
  {
    int size = 128 * 192;
    for (int channel = 0; channel < size; channel++)
    {
      std::string hname = (boost::format("h_cemc_channel_pedestal_%d") % channel).str();
      h_cemc_channel_pedestal[channel] = new TH1F(hname.c_str(), hname.c_str(), 2000, -0.5, 2000.5);
      h_cemc_channel_pedestal[channel]->SetDirectory(nullptr);

      std::string hnameE = (boost::format("h_cemc_channel_energy_%d") % channel).str();
      h_cemc_channel_energy[channel] = new TH1F(hnameE.c_str(), hnameE.c_str(), 1000, -50, 50);
      h_cemc_channel_energy[channel]->SetDirectory(nullptr);
    }
  }
  //--------OHCal--------//
  {
    int size = 32 * 48;
    for (int channel = 0; channel < size; channel++)
    {
      std::string hname = (boost::format("h_ohcal_channel_pedestal_%d") % channel).str();
      h_ohcal_channel_pedestal[channel] = new TH1F(hname.c_str(), hname.c_str(), 2000, -0.5, 2000.5);
      h_ohcal_channel_pedestal[channel]->SetDirectory(nullptr);

      std::string hnameE = (boost::format("h_ohcal_channel_energy_%d") % channel).str();
      h_ohcal_channel_energy[channel] = new TH1F(hnameE.c_str(), hnameE.c_str(), 1000, -50, 50);
      h_ohcal_channel_energy[channel]->SetDirectory(nullptr);
    }
  }
  //--------IHCal-------//
  {
    int size = 32 * 48;
    for (int channel = 0; channel < size; channel++)
    {
      std::string hname = (boost::format("h_ihcal_channel_pedestal_%d") % channel).str();
      h_ihcal_channel_pedestal[channel] = new TH1F(hname.c_str(), hname.c_str(), 2000, -0.5, 2000.5);
      h_ihcal_channel_pedestal[channel]->SetDirectory(nullptr);

      std::string hnameE = (boost::format("h_ihcal_channel_energy_%d") % channel).str();
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
  // h_triggerVec = new TH1F("h_CaloValid_triggerVec", "", 64, 0, 64); commented out due to memory allocation issue w/ redef from line 1009
  pr_ldClus_trig =
      new TProfile("pr_CaloValid_ldClus_trig", "", 64, 0, 64, 0, 10);
  for (int i = 0; i < 64; i++)
  {
    if (!(std::find(trigOfInterest.begin(), trigOfInterest.end(), i) != trigOfInterest.end()))
    {
      continue;
    }
    h_edist[i] = new TH2F(
        boost::str(boost::format("h_CaloValid_edist_trig%d") % i).c_str(), "",
        64, -1.2, 1.2, 128, -3.1415, 3.1415);
    h_ldClus_trig[i] = new TH1F(
        boost::str(boost::format("h_CaloValid_ldClus_trig%d") % i).c_str(), "",
        18, 1, 10);
    pr_evtNum_ldClus_trig[i] = new TProfile(
        boost::str(boost::format("pr_CaloValid_evtNum_ldClus_trig%d") % i)
            .c_str(),
        "", 100000, 0, 100000, 0, 10);
    pr_rejection[i] = new TProfile(
        boost::str(boost::format("pr_CaloValid_rejection_trig%d") % i).c_str(),
        "", 100000, 0, 100000, 0, 50000);
    pr_livetime[i] = new TProfile(
        boost::str(boost::format("pr_CaloValid_livetime_trig%d") % i).c_str(),
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

