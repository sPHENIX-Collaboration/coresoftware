#include "CaloValid.h"

// Calo includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <qautils/QAHistManagerDef.h>

#include <boost/format.hpp>

// MBD
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include <TTree.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <cmath>     // for log10, pow, sqrt, abs, M_PI
#include <iostream>  // for operator<<, endl, basic_...
#include <limits>
#include <map>  // for operator!=, _Rb_tree_con...
#include <string>
#include <utility>  // for pair

CaloValid::CaloValid(const std::string& name)
  : SubsysReco(name)
  , detector("HCALIN")
{
}

CaloValid::~CaloValid() = default;

int CaloValid::Init(PHCompositeNode* /*unused*/)
{
  if (m_debug)
  {
    std::cout << "In CaloValid::Init" << std::endl;
  }

  createHistos();

  if (m_debug)
  {
    std::cout << "Leaving CaloValid::Init" << std::endl;
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
  if (m_debug)
  {
    std::cout << _eventcounter << std::endl;
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  float totalcemc = 0.;
  float totalihcal = 0.;
  float totalohcal = 0.;
  float totalmbd = 0.;
  float totalzdc = 0.;
  float totalzdcsouthraw = 0.;
  float totalzdcnorthraw = 0.;
  float totalzdcsouthcalib = 0.;
  float totalzdcnorthcalib = 0.;

  float emcaldownscale = 1000000. / 800.;
  float ihcaldownscale = 40000. / 300.;
  float ohcaldownscale = 250000. / 600.;
  float mbddownscale = 2000.0;
  float zdcdownscale = 1e4;

  float adc_threshold = 15.;

  float emcal_hit_threshold = 0.5;  // GeV
  float ohcal_hit_threshold = 0.5;
  float ihcal_hit_threshold = 0.25;

  int max_zdc_t = -1;
  int max_emcal_t = -1;
  int max_ihcal_t = -1;
  int max_ohcal_t = -1;

  // get time estimate
  auto hzdctime_cut = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%szdctime_cut") % getHistoPrefix()).c_str()));
  auto hemcaltime_cut = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%semcaltime_cut") % getHistoPrefix()).c_str()));
  auto hihcaltime_cut = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sihcaltime_cut") % getHistoPrefix()).c_str()));
  auto hohcaltime_cut = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sohcaltime_cut") % getHistoPrefix()).c_str()));

  max_zdc_t = Getpeaktime(hzdctime_cut);
  max_emcal_t = Getpeaktime(hemcaltime_cut);
  max_ihcal_t = Getpeaktime(hihcaltime_cut);
  max_ohcal_t = Getpeaktime(hohcaltime_cut);

  //----------------------------------vertex------------------------------------------------------//
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << "CaloValid GlobalVertexMap node is missing" << std::endl;
  }
  float vtx_z = std::numeric_limits<float>::quiet_NaN();
  auto hvtx_z_raw = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%svtx_z_raw") % getHistoPrefix()).c_str()));

  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex* vtx = vertexmap->begin()->second;
    if (vtx)
    {
      vtx_z = vtx->get_z();
    }
    hvtx_z_raw->Fill(vtx_z);
  }

  //---------------------------calibrated towers-------------------------------//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      auto h_cemc_e_chi2 = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%scemc_e_chi2") % getHistoPrefix()).c_str()));
      auto h_cemc_status = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%scemc_status") % getHistoPrefix()).c_str()));
      auto hemcaltime = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%semcaltime") % getHistoPrefix()).c_str()));
      auto h_cemc_etaphi_time = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%scemc_etaphi_time") % getHistoPrefix()).c_str()));
      auto h_cemc_etaphi = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%scemc_etaphi") % getHistoPrefix()).c_str()));
      auto h_cemc_etaphi_wQA = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%scemc_etaphi_wQA") % getHistoPrefix()).c_str()));
      auto h_cemc_etaphi_badChi2 = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%scemc_etaphi_badChi2") % getHistoPrefix()).c_str()));

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
        hemcaltime_cut->Fill(_time);
        bool isGood = tower->get_isGood();
        uint8_t status = tower->get_status();
        for (int is = 0; is < 8; is++)
        {
          if (status & 1U)  // clang-tidy mark 1 as unsigned
          {
            h_cemc_status->Fill(is);
          }
          status = status >> 1U;  // clang-tidy mark 1 as unsigned
        }
        if (_time > (max_emcal_t - _range) && _time < (max_emcal_t + _range))
        {
          totalcemc += offlineenergy;
          hemcaltime->Fill(_time);
          if (offlineenergy > emcal_hit_threshold)
          {
            h_cemc_etaphi_time->Fill(ieta, iphi, _timef);
            h_cemc_etaphi->Fill(ieta, iphi);
            if (isGood)
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
        }
      }
    }
  }

  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    if (towers)
    {
      auto h_ihcal_e_chi2 = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sihcal_e_chi2") % getHistoPrefix()).c_str()));
      auto h_ihcal_status = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sihcal_status") % getHistoPrefix()).c_str()));
      auto hihcaltime = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sihcaltime") % getHistoPrefix()).c_str()));
      auto h_hcalin_etaphi_time = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sihcal_etaphi_time") % getHistoPrefix()).c_str()));
      auto h_hcalin_etaphi = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sihcal_etaphi") % getHistoPrefix()).c_str()));
      auto h_hcalin_etaphi_wQA = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sihcal_etaphi_wQA") % getHistoPrefix()).c_str()));
      auto h_hcalin_etaphi_badChi2 = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sihcal_etaphi_badChi2") % getHistoPrefix()).c_str()));

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
        hihcaltime_cut->Fill(_time);
        h_ihcal_e_chi2->Fill(offlineenergy, tower->get_chi2());
        bool isGood = tower->get_isGood();
        h_ihcal_status->Fill(tower->get_status());

        uint8_t status = tower->get_status();
        for (int is = 0; is < 8; is++)
        {
          if (status & 1U)  // clang-tidy mark 1 as unsigned
          {
            h_ihcal_status->Fill(is);
          }
          status = status >> 1U;  // clang-tidy mark 1 as unsigned
        }

        if (_time > (max_ihcal_t - _range) && _time < (max_ihcal_t + _range))
        {
          totalihcal += offlineenergy;
          hihcaltime->Fill(_time);

          if (offlineenergy > ihcal_hit_threshold)
          {
            h_hcalin_etaphi->Fill(ieta, iphi);
            h_hcalin_etaphi_time->Fill(ieta, iphi, _timef);
            if (isGood)
            {
              h_hcalin_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
            }
            if (tower->get_isBadChi2())
            {
              h_hcalin_etaphi_badChi2->Fill(ieta, iphi, 1);
            }
            else
            {
              h_hcalin_etaphi_badChi2->Fill(ieta, iphi, 0);
            }
          }
        }
      }
    }
  }

  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    if (towers)
    {
      auto h_ohcal_e_chi2 = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sohcal_e_chi2") % getHistoPrefix()).c_str()));
      auto h_ohcal_status = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sohcal_status") % getHistoPrefix()).c_str()));
      auto hohcaltime = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sohcaltime") % getHistoPrefix()).c_str()));
      auto h_hcalout_etaphi_time = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sohcal_etaphi_time") % getHistoPrefix()).c_str()));
      auto h_hcalout_etaphi = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sohcal_etaphi") % getHistoPrefix()).c_str()));
      auto h_hcalout_etaphi_wQA = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sohcal_etaphi_wQA") % getHistoPrefix()).c_str()));
      auto h_hcalout_etaphi_badChi2 = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sohcal_etaphi_badChi2") % getHistoPrefix()).c_str()));

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
        hohcaltime_cut->Fill(_time);
        h_ohcal_e_chi2->Fill(offlineenergy, tower->get_chi2());
        bool isGood = tower->get_isGood();
        h_ohcal_status->Fill(tower->get_status());

        uint8_t status = tower->get_status();
        for (int is = 0; is < 8; is++)
        {
          if (status & 1U)  // clang-tidy mark 1 as unsigned
          {
            h_ohcal_status->Fill(is);
          }
          status = status >> 1U;  // clang-tidy mark 1 as unsigned
        }

        if (_time > (max_ohcal_t - _range) && _time < (max_ohcal_t + _range))
        {
          totalohcal += offlineenergy;
          hohcaltime->Fill(_time);

          if (offlineenergy > ohcal_hit_threshold)
          {
            h_hcalout_etaphi->Fill(ieta, iphi);
            h_hcalout_etaphi_time->Fill(ieta, iphi, _timef);
            if (isGood)
            {
              h_hcalout_etaphi_wQA->Fill(ieta, iphi, offlineenergy);
            }
            if (tower->get_isBadChi2())
            {
              h_hcalout_etaphi_badChi2->Fill(ieta, iphi, 1);
            }
            else
            {
              h_hcalout_etaphi_badChi2->Fill(ieta, iphi, 0);
            }
          }
        }
      }
    }
  }

  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_ZDC");
    if (towers)
    {
      auto hzdctime = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%szdctime") % getHistoPrefix()).c_str()));

      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        float offlineenergy = tower->get_energy();
        int _time = towers->get_tower_at_channel(channel)->get_time();
        hzdctime_cut->Fill(_time);
        if (channel == 0 || channel == 2 || channel == 4)
        {
          totalzdcsouthcalib += offlineenergy;
        }
        if (channel == 8 || channel == 10 || channel == 12)
        {
          totalzdcnorthcalib += offlineenergy;
        }
        if (channel == 0 || channel == 2 || channel == 4 || channel == 8 || channel == 10 || channel == 12)
        {
          if (_time > (max_zdc_t - _range) && _time < (max_zdc_t + _range))
          {
            totalzdc += offlineenergy;
            hzdctime->Fill(_time);
          }
        }
      }
    }
  }

  //-------------------------- raw tower ------------------------------//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        float offlineenergy = tower->get_energy();
        if (channel == 0 || channel == 2 || channel == 4)
        {
          totalzdcsouthraw += offlineenergy;
        }
        if (channel == 8 || channel == 10 || channel == 12)
        {
          totalzdcnorthraw += offlineenergy;
        }
      }
    }
  }

  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
    if (towers)
    {
      auto h_cemc_etaphi_fracHitADC = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%scemc_etaphi_fracHitADC") % getHistoPrefix()).c_str()));

      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);

        float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_cemc_etaphi_fracHitADC->Fill(ieta, iphi, 1);
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
      auto h_hcalout_etaphi_fracHitADC = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sohcal_etaphi_fracHitADC") % getHistoPrefix()).c_str()));

      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);

        float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_hcalout_etaphi_fracHitADC->Fill(ieta, iphi, 1);
        }
        else
        {
          h_hcalout_etaphi_fracHitADC->Fill(ieta, iphi, 0);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
    if (towers)
    {
      auto h_hcalin_etaphi_fracHitADC = dynamic_cast<TProfile2D*>(hm->getHisto(boost::str(boost::format("%sihcal_etaphi_fracHitADC") % getHistoPrefix()).c_str()));

      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);

        float raw_energy = tower->get_energy();
        if (raw_energy > adc_threshold)
        {
          h_hcalin_etaphi_fracHitADC->Fill(ieta, iphi, 1);
        }
        else
        {
          h_hcalin_etaphi_fracHitADC->Fill(ieta, iphi, 0);
        }
      }
    }
  }

  //--------------------------- MBD ----------------------------------------//
  MbdPmtContainer* bbcpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!bbcpmts)
  {
    std::cout << "makeMBDTrees::process_event: Could not find MbdPmtContainer, aborting" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  int nPMTs = bbcpmts->get_npmt();
  for (int i = 0; i < nPMTs; i++)
  {
    MbdPmtHit* mbdpmt = bbcpmts->get_pmt(i);
    float pmtadc = mbdpmt->get_q();
    totalmbd += pmtadc;
  }
  auto h_emcal_mbd_correlation = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%semcal_mbd_correlation") % getHistoPrefix()).c_str()));
  auto h_ihcal_mbd_correlation = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sihcal_mbd_correlation") % getHistoPrefix()).c_str()));
  auto h_ohcal_mbd_correlation = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%sohcal_mbd_correlation") % getHistoPrefix()).c_str()));
  auto h_emcal_hcal_correlation = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%semcal_hcal_correlation") % getHistoPrefix()).c_str()));
  auto h_emcal_zdc_correlation = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%szdc_emcal_correlation") % getHistoPrefix()).c_str()));
  auto h_totalzdc_e = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%stotalzdc_e") % getHistoPrefix()).c_str()));
  auto hzdcSouthraw = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%szdcSouthraw") % getHistoPrefix()).c_str()));
  auto hzdcNorthraw = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%szdcNorthraw") % getHistoPrefix()).c_str()));
  auto hzdcSouthcalib = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%szdcSouthcalib") % getHistoPrefix()).c_str()));
  auto hzdcNorthcalib = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%szdcNorthcalib") % getHistoPrefix()).c_str()));

  h_emcal_mbd_correlation->Fill(totalcemc / emcaldownscale, totalmbd / mbddownscale);
  h_ihcal_mbd_correlation->Fill(totalihcal / ihcaldownscale, totalmbd / mbddownscale);
  h_ohcal_mbd_correlation->Fill(totalohcal / ohcaldownscale, totalmbd / mbddownscale);
  h_emcal_hcal_correlation->Fill(totalcemc / emcaldownscale, totalohcal / ohcaldownscale);
  h_emcal_zdc_correlation->Fill(totalcemc / emcaldownscale, totalzdc / zdcdownscale);
  h_totalzdc_e->Fill(totalzdc);

  hzdcSouthraw->Fill(totalzdcsouthraw);
  hzdcNorthraw->Fill(totalzdcnorthraw);
  hzdcSouthcalib->Fill(totalzdcsouthcalib);
  hzdcNorthcalib->Fill(totalzdcnorthcalib);

  //------------------------------ clusters & pi0 ------------------------------//
  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_POS_COR_CEMC");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "funkyCaloStuff::process_event - Fatal Error - CLUSTER_CEMC node is missing. " << std::endl;
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
  auto h_clusE = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sclusE") % getHistoPrefix()).c_str()));
  auto h_etaphi_clus = dynamic_cast<TH2*>(hm->getHisto(boost::str(boost::format("%setaphi_clus") % getHistoPrefix()).c_str()));
  auto h_InvMass = dynamic_cast<TH1*>(hm->getHisto(boost::str(boost::format("%sInvMass") % getHistoPrefix()).c_str()));

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
        h_InvMass->Fill(pi0.M());
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::End(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloValid::Getpeaktime(TH1* h)
{
  int getmaxtime, tcut = -1;

  for (int bin = 1; bin < h->GetNbinsX() + 1; bin++)
  {
    double c = h->GetBinContent(bin);
    double max = h->GetMaximum();
    int bincenter = h->GetBinCenter(bin);
    if (max == c)
    {
      getmaxtime = bincenter;
      if (getmaxtime != -1)
      {
        tcut = getmaxtime;
      }
    }
  }

  return tcut;
}

TH2F* CaloValid::LogYHist2D(const std::string& name, const std::string& title, int xbins_in, double xmin, double xmax, int ybins_in, double ymin, double ymax)
{
  Double_t logymin = std::log10(ymin);
  Double_t logymax = std::log10(ymax);
  Double_t binwidth = (logymax - logymin) / ybins_in;
  Double_t ybins[ybins_in + 1];

  for (Int_t i = 0; i <= ybins_in + 1; i++)
  {
    ybins[i] = pow(10, logymin + i * binwidth);
  }

  TH2F* h = new TH2F(name.c_str(), title.c_str(), xbins_in, xmin, xmax, ybins_in, ybins);

  return h;
}
std::string CaloValid::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void CaloValid::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create and register your histos (all types) here
  {
    auto h = new TH2F(boost::str(boost::format("%semcal_mbd_correlation") % getHistoPrefix()).c_str(), ";emcal;mbd", 100, 0, 1, 100, 0, 1);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%sohcal_mbd_correlation") % getHistoPrefix()).c_str(), ";ohcal;mbd", 100, 0, 1, 100, 0, 1);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%sihcal_mbd_correlation") % getHistoPrefix()).c_str(), ";ihcal;mbd", 100, 0, 1, 100, 0, 1);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%semcal_hcal_correlation") % getHistoPrefix()).c_str(), ";emcal;hcal", 100, 0, 1, 100, 0, 1);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%scemc_etaphi") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%sihcal_etaphi") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%sohcal_etaphi") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%scemc_etaphi_wQA") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%sihcal_etaphi_wQA") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%sohcal_etaphi_wQA") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(boost::str(boost::format("%sihcal_status") % getHistoPrefix()).c_str(), "", 256, 0, 256);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(boost::str(boost::format("%sohcal_status") % getHistoPrefix()).c_str(), "", 256, 0, 256);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(boost::str(boost::format("%scemc_status") % getHistoPrefix()).c_str(), "", 256, 0, 256);
    hm->registerHisto(h);
  }

  {
    auto h = LogYHist2D(boost::str(boost::format("%scemc_e_chi2") % getHistoPrefix()).c_str(), "", 500, -2, 30, 1000, 0.5, 5e6);
    hm->registerHisto(h);
  }
  {
    auto h = LogYHist2D(boost::str(boost::format("%sihcal_e_chi2") % getHistoPrefix()).c_str(), "", 500, -2, 30, 1000, 0.5, 5e6);
    hm->registerHisto(h);
  }
  {
    auto h = LogYHist2D(boost::str(boost::format("%sohcal_e_chi2") % getHistoPrefix()).c_str(), "", 500, -2, 30, 1000, 0.5, 5e6);
    hm->registerHisto(h);
  }

  {
    auto h = new TProfile2D(boost::str(boost::format("%scemc_etaphi_time") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_time") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_time") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%scemc_etaphi_fracHitADC") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_fracHitADC") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_fracHitADC") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%scemc_etaphi_badChi2") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_badChi2") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
    hm->registerHisto(h);
  }
  {
    auto h = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_badChi2") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
    hm->registerHisto(h);
  }
  // 1D distributions
  {
    auto h = new TH1F(boost::str(boost::format("%sInvMass") % getHistoPrefix()).c_str(), "Invariant Mass", 120, 0, 1.2);
    hm->registerHisto(h);
  }
  // ZDC QA plots
  {
    auto h = new TH1D(boost::str(boost::format("%szdcSouthraw") % getHistoPrefix()).c_str(), "hzdcSouthraw", 1500, 0, 15000);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%szdcNorthraw") % getHistoPrefix()).c_str(), "hzdcNorthraw", 1500, 0, 15000);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%szdcSouthcalib") % getHistoPrefix()).c_str(), "hzdcSouthcalib", 1500, 0, 15000);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%szdcNorthcalib") % getHistoPrefix()).c_str(), "hzdcNorthcalib", 1500, 0, 15000);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%stotalzdc_e") % getHistoPrefix()).c_str(), "", 200, 0, 2e4);
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(boost::str(boost::format("%szdc_emcal_correlation") % getHistoPrefix()).c_str(), ";emcal;zdc", 100, 0, 1, 100, 0, 1);
    hm->registerHisto(h);
  }
  // vertex distributions
  {
    auto h = new TH1D(boost::str(boost::format("%svtx_z_raw") % getHistoPrefix()).c_str(), "hvtx_z_raw", 201, -100.5, 100.5);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%svtx_z_cut") % getHistoPrefix()).c_str(), "hvtx_z_cut", 201, -100.5, 100.5);
    hm->registerHisto(h);
  }

  // raw timing information
  {
    auto h = new TH1D(boost::str(boost::format("%szdctime_cut") % getHistoPrefix()).c_str(), "hzdctime_cut", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%semcaltime_cut") % getHistoPrefix()).c_str(), "hemcaltime_cut", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%sihcaltime_cut") % getHistoPrefix()).c_str(), "hihcaltime_cut", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%sohcaltime_cut") % getHistoPrefix()).c_str(), "hohcaltime_cut", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }

  // extracted timing information
  {
    auto h = new TH1D(boost::str(boost::format("%szdctime") % getHistoPrefix()).c_str(), "hzdctime", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%semcaltime") % getHistoPrefix()).c_str(), "hemcaltime", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%sihcaltime") % getHistoPrefix()).c_str(), "hihcaltime", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1D(boost::str(boost::format("%sohcaltime") % getHistoPrefix()).c_str(), "hohcaltime", 50, -17.5, 32.5);
    hm->registerHisto(h);
  }

  // cluster QA
  {
    auto h = new TH2F(boost::str(boost::format("%setaphi_clus") % getHistoPrefix()).c_str(), "", 140, -1.2, 1.2, 64, -1 * M_PI, M_PI);
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(boost::str(boost::format("%sclusE") % getHistoPrefix()).c_str(), "", 100, 0, 10);
    hm->registerHisto(h);
  }
}