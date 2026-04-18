///////////////////////
//EMCal Shower Shape QA
//
///////////////////////
#include "EMCalShowerShapes.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <calotrigger/TriggerAnalyzer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <qautils/QAHistManagerDef.h>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

namespace
{
  void shift_tower_index(int& ieta, int& iphi, int etadiv, int phidiv)
  {
    while (iphi < 0)
    {
      iphi += phidiv;
    }
    while (iphi >= phidiv)
    {
      iphi -= phidiv;
    }
    if (ieta < 0 || ieta >= etadiv)
    {
      ieta = -1;
    }
  }
}

EMCalShowerShapes::EMCalShowerShapes(const std::string &modulename, const std::string &inputnode, const std::string &histtag)
  : SubsysReco(modulename)
  , m_modulename(modulename)
  , m_inputnode(inputnode)
  , m_histtag(histtag)
  , m_trgToSelect(JetQADefs::GL1::MBDNSPhoton1)
  , m_doTrgSelect(false)
{
}

EMCalShowerShapes::~EMCalShowerShapes()
{
  delete m_analyzer;
}

int EMCalShowerShapes::Init(PHCompositeNode* /*topNode*/)
{
  delete m_analyzer;
  m_analyzer = new TriggerAnalyzer();

  m_manager = QAHistManagerDef::getHistoManager();
  if (!m_manager)
  {
    std::cerr << PHWHERE << "PANIC: couldn't grab histogram manager!" << std::endl;
    gSystem->Exit(1);
  }

  std::string smallModuleName = m_modulename;
  std::transform(smallModuleName.begin(), smallModuleName.end(), smallModuleName.begin(), ::tolower);

  std::vector<std::string> vecHistNames = {
      "cluster_et",
      "e11_to_e33",
      "e33_to_e55",
      "e55_to_e77",
      "e32_to_e35",
      "weta",
      "wphi",
      "weta_cogx",
      "wphi_cogx",
      "detamax",
      "dphimax",
      "mean_time",
      "iso04_emcal",
      "weta_vs_et",
      "wphi_vs_et"};

  for (auto &histName : vecHistNames)
  {
    histName.insert(0, "h_" + smallModuleName + "_");
    if (!m_histtag.empty())
    {
      histName.append("_" + m_histtag);
    }
  }

  h_cluster_et = new TH1F(vecHistNames[0].data(), "", 120, 0, 30);
  h_cluster_et->GetXaxis()->SetTitle("E_{T} [GeV]");

  h_e11oe33 = new TH1F(vecHistNames[1].data(), "", 26, -0.02, 1.02);
  h_e11oe33->GetXaxis()->SetTitle("e11/e33");

  h_e33oe55 = new TH1F(vecHistNames[2].data(), "", 26, -0.02, 1.02);
  h_e33oe55->GetXaxis()->SetTitle("e33/e55");

  h_e55oe77 = new TH1F(vecHistNames[3].data(), "", 26, -0.02, 1.02);
  h_e55oe77->GetXaxis()->SetTitle("e55/e77");
  
  h_e32oe35 = new TH1F(vecHistNames[4].data(), "", 26, -0.02, 1.02);
  h_e32oe35->GetXaxis()->SetTitle("e32/e35");

  h_weta = new TH1F(vecHistNames[5].data(), "", 120, 0, 2);
  h_weta->GetXaxis()->SetTitle("w#eta");

  h_wphi = new TH1F(vecHistNames[6].data(), "", 120, 0, 2);
  h_wphi->GetXaxis()->SetTitle("w#phi");

  h_weta_cogx = new TH1F(vecHistNames[7].data(), "", 50, 0, 2);
  h_weta_cogx->GetXaxis()->SetTitle("w#eta_cogx");

  h_wphi_cogx = new TH1F(vecHistNames[8].data(), "", 50, 0, 2);
  h_wphi_cogx->GetXaxis()->SetTitle("w#phi_cogx");

  h_detamax = new TH1F(vecHistNames[9].data(), "", 10, -0.5, 9.5);
  h_detamax->GetXaxis()->SetTitle("detamax");

  h_dphimax = new TH1F(vecHistNames[10].data(), "", 20, -0.5, 19.5);
  h_dphimax->GetXaxis()->SetTitle("dphimax");

  h_mean_time = new TH1F(vecHistNames[11].data(), "", 200, -20, 20);
  h_mean_time->GetXaxis()->SetTitle("cluster mean time");

  h_iso04_emcal = new TH1F(vecHistNames[12].data(), "", 200, -10, 40);
  h_iso04_emcal->GetXaxis()->SetTitle("iso_{0.4}^{EMCal} [GeV]");

  h_weta_vs_et = new TH2F(vecHistNames[13].data(), "", 120, 0, 30, 120, 0, 6);
  h_weta_vs_et->GetXaxis()->SetTitle("E_{T} [GeV]");
  h_weta_vs_et->GetYaxis()->SetTitle("w#eta");

  h_wphi_vs_et = new TH2F(vecHistNames[14].data(), "", 120, 0, 30, 120, 0, 6);
  h_wphi_vs_et->GetXaxis()->SetTitle("E_{T} [GeV]");
  h_wphi_vs_et->GetYaxis()->SetTitle("w#phi");

  // Register histograms here to preserve them even if files are closedß
  m_manager->registerHisto(h_cluster_et);
  m_manager->registerHisto(h_e11oe33);
  m_manager->registerHisto(h_e33oe55);
  m_manager->registerHisto(h_e55oe77);
  m_manager->registerHisto(h_e32oe35);
  m_manager->registerHisto(h_weta);
  m_manager->registerHisto(h_wphi);
  m_manager->registerHisto(h_weta_cogx);
  m_manager->registerHisto(h_wphi_cogx);
  m_manager->registerHisto(h_detamax);
  m_manager->registerHisto(h_dphimax);
  m_manager->registerHisto(h_mean_time);
  m_manager->registerHisto(h_iso04_emcal);
  m_manager->registerHisto(h_weta_vs_et);
  m_manager->registerHisto(h_wphi_vs_et);

  return Fun4AllReturnCodes::EVENT_OK;
}

int EMCalShowerShapes::InitRun(PHCompositeNode* topNode)
{
  if (!LoadEMCalNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

bool EMCalShowerShapes::LoadEMCalNodes(PHCompositeNode *topNode)
{
  m_emc_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  m_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");

  const bool have_nodes = (m_emc_tower_container && m_geomEM);
  if (!have_nodes && !m_reportedMissingCaloNodes)
  {
    std::cout << PHWHERE << "EMCalShowerShapes::LoadEMCalNodes - missing TOWERINFO_CALIB_CEMC or TOWERGEOM_CEMC" << std::endl;
    m_reportedMissingCaloNodes = true;
  }
  return have_nodes;
}

float EMCalShowerShapes::GetVertexZ(PHCompositeNode *topNode) const
{
  MbdVertexMap* vertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  if (!vertexmap || vertexmap->empty())
  {
    return 0.0F;
  }

  MbdVertex* vtx = vertexmap->begin()->second;
  if (!vtx)
  {
    return 0.0F;
  }

  return vtx->get_z();
}

int EMCalShowerShapes::process_event(PHCompositeNode *topNode)
{
  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, m_inputnode);
  if (!clusterContainer)
  {
    if (!m_reportedMissingClusterNode)
    {
      std::cout << PHWHERE << "EMCalShowerShapes::process_event - missing node " << m_inputnode << std::endl;
      m_reportedMissingClusterNode = true;
    }
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!LoadEMCalNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_doTrgSelect)
  {
    m_analyzer->decodeTriggers(topNode);
    if (!JetQADefs::DidTriggerFire(m_trgToSelect, m_analyzer))
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  const float vertex_z = GetVertexZ(topNode);
  if (m_doMbdZvtxCut && std::abs(vertex_z) >= m_mbdZvtxMax)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const CLHEP::Hep3Vector vertex_vec(0, 0, vertex_z);

  RawClusterContainer::ConstRange clusters = clusterContainer->getClusters();
  for (auto iter = clusters.first; iter != clusters.second; ++iter)
  {
    RawCluster* cluster = iter->second;
    if (!cluster)
    {
      continue;
    }

    const float eta = RawClusterUtility::GetPseudorapidity(*cluster, vertex_vec);
    if (m_doClusterEtaCut && std::abs(eta) >= m_clusterEtaMax)
    {
      continue;
    }

    const float phi = RawClusterUtility::GetAzimuthAngle(*cluster, vertex_vec);
    const float et = cluster->get_energy() / std::cosh(eta);
    if (m_doClusterETCut && et < m_clusterETMin)
    {
      continue;
    }

    ShowerShapeData data;
    if (!CalculateShowerShapes(cluster, eta, phi, et, vertex_z, data))
    {
      continue;
    }

    h_cluster_et->Fill(et);
    if (data.e33 > 0) {
      h_e11oe33->Fill(data.e11 / data.e33);
    }
    if (data.e55 > 0) {
      h_e33oe55->Fill(data.e33 / data.e55);
    }
    if (data.e77 > 0) {
      h_e55oe77->Fill(data.e55 / data.e77);
    }
    if (data.e35 > 0) {
      h_e32oe35->Fill(data.e32 / data.e35);
    }
    h_weta->Fill(data.weta);
    h_wphi->Fill(data.wphi);
    h_weta_cogx->Fill(data.weta_cogx);
    h_wphi_cogx->Fill(data.wphi_cogx);
    h_detamax->Fill(data.detamax);
    h_dphimax->Fill(data.dphimax);
    h_mean_time->Fill(data.mean_time);
    h_iso04_emcal->Fill(data.iso04_emcal);
    h_weta_vs_et->Fill(et, data.weta);
    h_wphi_vs_et->Fill(et, data.wphi);

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool EMCalShowerShapes::CalculateShowerShapes(RawCluster* cluster, float cluster_eta, float cluster_phi, float cluster_et, float vertex_z, ShowerShapeData& data) const
{
  std::vector<float> showershape = cluster->get_shower_shapes(m_shape_min_tower_E);
  if (showershape.empty())
  {
    return false;
  }

  const std::pair<int, int> leadtowerindex = cluster->get_lead_tower();
  const int lead_ieta = leadtowerindex.first;
  const int lead_iphi = leadtowerindex.second;

  const float avg_eta = showershape[4] + 0.5F;
  const float avg_phi = showershape[5] + 0.5F;
  const int maxieta = std::floor(avg_eta);
  const int maxiphi = std::floor(avg_phi);

  int detamax = 0;
  int dphimax = 0;
  float clusteravgtime = 0.0F;
  float cluster_total_e = 0.0F;
  const RawCluster::TowerMap& tower_map = cluster->get_towermap();
  std::set<unsigned int> towers_in_cluster;
  for (auto tower_iter : tower_map)
  {
    RawTowerDefs::keytype tower_key = tower_iter.first;
    const int ieta = RawTowerDefs::decode_index1(tower_key);
    const int iphi = RawTowerDefs::decode_index2(tower_key);

    const unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, iphi);
    towers_in_cluster.insert(towerinfokey);
    TowerInfo* towerinfo = m_emc_tower_container->get_tower_at_key(towerinfokey);
    if (towerinfo)
    {
      clusteravgtime += towerinfo->get_time() * towerinfo->get_energy();
      cluster_total_e += towerinfo->get_energy();
    }

    int totalphibins = 256;
    auto dphiwrap = [totalphibins](int towerphi, int maxiphi_arg)
    {
      int idphi = towerphi - maxiphi_arg;
      if (idphi > totalphibins / 2)
      {
        idphi -= totalphibins;
      }
      if (idphi < -totalphibins / 2)
      {
        idphi += totalphibins;
      }
      return idphi;
    };

    const int deta = ieta - lead_ieta;
    const int dphi_val = dphiwrap(iphi, lead_iphi);
    detamax = std::max(std::abs(deta), detamax);
    dphimax = std::max(std::abs(dphi_val), dphimax);
  }

  if (cluster_total_e > 0)
  {
    clusteravgtime /= cluster_total_e;
  }
  else
  {
    std::cout << "cluster_total_e is 0(this should not happen!!!), setting clusteravgtime to NaN" << std::endl;
    clusteravgtime = std::numeric_limits<float>::quiet_NaN();
  }

  float E77[7][7] = {{0.0F}};
  int E77_ownership[7][7] = {{0}};

  for (int ieta = maxieta - 3; ieta < maxieta + 4; ++ieta)
  {
    for (int iphi = maxiphi - 3; iphi < maxiphi + 4; ++iphi)
    {
      if (ieta < 0 || ieta > 95)
      {
        E77[ieta - maxieta + 3][iphi - maxiphi + 3] = 0.0F;
        E77_ownership[ieta - maxieta + 3][iphi - maxiphi + 3] = 0;
        continue;
      }

      int temp_ieta = ieta;
      int temp_iphi = iphi;
      shift_tower_index(temp_ieta, temp_iphi, 96, 256);
      if (temp_ieta < 0)
      {
        continue;
      }

      const unsigned int towerinfokey = TowerInfoDefs::encode_emcal(temp_ieta, temp_iphi);
      //if (towers_in_cluster.find(towerinfokey) != towers_in_cluster.end())
      if (towers_in_cluster.contains(towerinfokey))
      {
        E77_ownership[ieta - maxieta + 3][iphi - maxiphi + 3] = 1;
      }

      TowerInfo* towerinfo = m_emc_tower_container->get_tower_at_key(towerinfokey);
      if (towerinfo && towerinfo->get_isGood())
      {
        const float energy = towerinfo->get_energy();
        if (energy > m_shape_min_tower_E)
        {
          E77[ieta - maxieta + 3][iphi - maxiphi + 3] = energy;
        }
      }
    }
  }

  float e11 = E77[3][3];
  float e32 = 0.0F;
  float e33 = 0.0F;
  float e35 = 0.0F;
  float e55 = 0.0F;
  float e77 = 0.0F;
  float weta = 0.0F;
  float wphi = 0.0F;
  float weta_cogx = 0.0F;
  float wphi_cogx = 0.0F;
  float Eetaphi = 0.0F;

  const float shift_eta = avg_eta - std::floor(avg_eta) - 0.5F;
  const float shift_phi = avg_phi - std::floor(avg_phi) - 0.5F;
  const float cog_eta = 3 + shift_eta;
  const float cog_phi = 3 + shift_phi;
  const int signphi = (avg_phi - std::floor(avg_phi)) > 0.5 ? 1 : -1;

  for (int i = 0; i < 7; ++i)
  {
    for (int j = 0; j < 7; ++j)
    {
      const int di = std::abs(i - 3);
      const int dj = std::abs(j - 3);
      const float di_float = i - cog_eta;
      const float dj_float = j - cog_phi;

      if (E77_ownership[i][j] == 1)
      {
        weta += E77[i][j] * di * di;
        wphi += E77[i][j] * dj * dj;
        Eetaphi += E77[i][j];
        if (i != 3 || j != 3)
        {
          weta_cogx += E77[i][j] * di_float * di_float;
          wphi_cogx += E77[i][j] * dj_float * dj_float;
        }
      }

      e77 += E77[i][j];
      if (di <= 1 && (dj == 0 || j == (3 + signphi)))
      {
        e32 += E77[i][j];
      }
      if (di <= 1 && dj <= 1)
      {
        e33 += E77[i][j];
      }
      if (di <= 1 && dj <= 2)
      {
        e35 += E77[i][j];
      }
      if (di <= 2 && dj <= 2)
      {
        e55 += E77[i][j];
      }
    }
  }

  if (Eetaphi > 0)
  {
    weta /= Eetaphi;
    wphi /= Eetaphi;
    weta_cogx /= Eetaphi;
    wphi_cogx /= Eetaphi;
  }
  /*else
  {
    weta = std::numeric_limits<float>::quiet_NaN();
    wphi = std::numeric_limits<float>::quiet_NaN();
    weta_cogx = std::numeric_limits<float>::quiet_NaN();
    wphi_cogx = std::numeric_limits<float>::quiet_NaN();
  }*/

  data.e11 = e11;
  data.e33 = e33;
  data.e32 = e32;
  data.e35 = e35;
  data.e55 = e55;
  data.e77 = e77;
  data.weta = weta;
  data.wphi = wphi;
  data.weta_cogx = weta_cogx;
  data.wphi_cogx = wphi_cogx;
  data.detamax = detamax;
  data.dphimax = dphimax;
  data.mean_time = clusteravgtime;
  data.iso04_emcal = CalculateLayerET(cluster_eta, cluster_phi, 0.4F, m_emc_tower_container, m_geomEM, vertex_z) - cluster_et;

  return true;
}

double EMCalShowerShapes::GetTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz) const
{
  if (!tower_geom)
  {
    return -9999;
  }
  if (vx == 0 && vy == 0 && vz == 0)
  {
    return tower_geom->get_eta();
  }

  const double radius = std::sqrt((tower_geom->get_center_x() - vx) * (tower_geom->get_center_x() - vx) +
                                  (tower_geom->get_center_y() - vy) * (tower_geom->get_center_y() - vy));
  const double theta = std::atan2(radius, tower_geom->get_center_z() - vz);
  return -std::log(std::tan(theta / 2.));
}

double EMCalShowerShapes::DeltaR(double eta1, double phi1, double eta2, double phi2) const
{
  double dphi = phi1 - phi2;
  while (dphi > M_PI)
  {
    dphi -= 2 * M_PI;
  }
  while (dphi <= -M_PI)
  {
    dphi += 2 * M_PI;
  }
  return std::sqrt(std::pow(eta1 - eta2, 2) + std::pow(dphi, 2));
}

float EMCalShowerShapes::CalculateLayerET(float seed_eta, float seed_phi, float radius, TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, float vertex_z) const
{
  if (!towerContainer || !geomContainer)
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float layer_et = 0.0F;
  const unsigned int ntowers = towerContainer->size();
  for (unsigned int channel = 0; channel < ntowers; ++channel)
  {
    TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
    if (!tower || !tower->get_isGood())
    {
      continue;
    }

    const unsigned int towerkey = towerContainer->encode_key(channel);
    const int ieta = towerContainer->getTowerEtaBin(towerkey);
    const int iphi = towerContainer->getTowerPhiBin(towerkey);

    const RawTowerDefs::keytype geom_key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
    RawTowerGeom* tower_geom = geomContainer->get_tower_geometry(geom_key);
    if (!tower_geom)
    {
      continue;
    }

    const double tower_eta = GetTowerEta(tower_geom, 0, 0, vertex_z);
    const double tower_phi = tower_geom->get_phi();
    if (DeltaR(seed_eta, seed_phi, tower_eta, tower_phi) >= radius)
    {
      continue;
    }

    const float energy = tower->get_energy();
    if (energy <= m_shape_min_tower_E)
    {
      continue;
    }

    layer_et += energy / std::cosh(tower_eta);
  }

  return layer_et;
}

//int EMCalShowerShapes::ResetEvent(PHCompositeNode* /*topNode*/)
//{
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//int EMCalShowerShapes::EndRun(const int /*runnumber*/)
//{
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//int EMCalShowerShapes::End(PHCompositeNode* /*topNode*/)
//{
//  return Fun4AllReturnCodes::EVENT_OK;
//}

//int EMCalShowerShapes::Reset(PHCompositeNode* /*topNode*/)
//{
//  return Fun4AllReturnCodes::EVENT_OK;
//}

void EMCalShowerShapes::Print(const std::string &what) const
{
  std::cout << "EMCalShowerShapes::Print(" << what << ")" << std::endl;
}
