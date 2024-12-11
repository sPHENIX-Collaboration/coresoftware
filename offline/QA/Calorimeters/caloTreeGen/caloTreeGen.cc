#include "caloTreeGen.h"

// For clusters and geometry
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>

// Tower stuff
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// GL1 Information
#include <ffarawobjects/Gl1Packet.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// ROOT stuff
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

// for cluster vertex correction
#include <CLHEP/Vector/ThreeVector.h>

#include <cstdint>
#include <iostream>
#include <map>
#include <utility>

//____________________________________________________________________________..
caloTreeGen::caloTreeGen(const std::string &name)
  : SubsysReco("CaloTreeGen")
  , T(nullptr)
  , Outfile(name)
{
  std::cout << "caloTreeGen::caloTreeGen(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
int caloTreeGen::Init(PHCompositeNode * /*topNode*/)
{
  delete out;  // make cppcheck happy (nullptrs can be deleted)
  out = new TFile(Outfile.c_str(), "RECREATE");

  T = new TTree("T", "T");

  // Electromagnetic Calorimeter
  if (storeEMCal)
  {
    T->Branch("emcTowE", &m_emcTowE);
    T->Branch("emcTowiEta", &m_emciEta);
    T->Branch("emcTowiPhi", &m_emciPhi);
    T->Branch("emcTime", &m_emcTime);
    T->Branch("emcChi2", &m_emcChi2);
    T->Branch("emcPed", &m_emcPed);

    // EMCal Cluster information
    if (storeEMCal && storeClusters)
    {
      T->Branch("clusterE", &m_clusterE);
      T->Branch("clusterPhi", &m_clusterPhi);
      T->Branch("clusterEta", &m_clusterEta);
      T->Branch("clusterPt", &m_clusterPt);
      T->Branch("clusterChi2", &m_clusterChi);
      T->Branch("clusterNtow", &m_clusterNtow);
      T->Branch("clusterTowMaxE", &m_clusterTowMaxE);
      T->Branch("clusterECore", &m_clusterECore);

      // Information for towers within clusters
      // Enabled by setting "DoFineClusters" in the macro
      if (storeEMCal && storeClusters && storeClusterDetails)
      {
        T->Branch("clusTowPhi", "vector<vector<int> >", &m_clusTowPhi);
        T->Branch("clusTowEta", "vector<vector<int> >", &m_clusTowEta);
        T->Branch("clusTowE", "vector<vector<float> >", &m_clusTowE);
      }
    }
  }
  // Outer Hadronic Calorimeter
  if (storeHCals)
  {
    T->Branch("ohcTowE", &m_ohcTowE);
    T->Branch("ohcTowiEta", &m_ohciTowEta);
    T->Branch("ohcTowiPhi", &m_ohciTowPhi);
    T->Branch("ohcTime", &m_ohcTime);
    T->Branch("ohcChi2", &m_ohcChi2);
    T->Branch("ohcPed", &m_ohcPed);

    // Inner Hadronic Calorimeter
    T->Branch("ihcTowE", &m_ihcTowE);
    T->Branch("ihcTowiEta", &m_ihciTowEta);
    T->Branch("ihcTowiPhi", &m_ihciTowPhi);
    T->Branch("ihcTime", &m_ihcTime);
    T->Branch("ihcChi2", &m_ihcChi2);
    T->Branch("ihcPed", &m_ihcPed);
  }
  // ZDC information
  if (storeZDC)
  {
    T->Branch("zdcTowE", &m_zdcTowE);
    T->Branch("zdcTowside", &m_zdcSide);

    // SMD information
    T->Branch("smdE", &m_smdE);
    T->Branch("smdSide", &m_smdSide);
  }
  // Total
  T->Branch("totalCaloEEMCal", &totalCaloEEMCal);
  T->Branch("totalCaloEOHCal", &totalCaloEOHCal);
  T->Branch("totalCaloEIHCal", &totalCaloEIHCal);
  T->Branch("totalCaloEZDC", &totalCaloEZDC);
  T->Branch("zvertex", &m_vertex);

  // GL1 Information
  if (storeTrig)
  {
    T->Branch("triggerVector", &m_triggerVector);
  }

  zVertex = new TH1F("zVertex", "zVertex", 200, -100, 100);

  // so that the histos actually get written out
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Print("NODETREE");

  std::cout << "caloTreeGen::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::process_event(PHCompositeNode *topNode)
{
  
  if(eventNum >= 100)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

  m_vertex = -9999;
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << "GlobalVertexMap node is missing" << std::endl;
  }
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx)
    {
      m_vertex = vtx->get_z();
      zVertex->Fill(m_vertex);
    }
  }

  // Information on clusters
  RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, m_clusterNode.c_str());
  if (!clusterContainer && storeClusters)
  {
    std::cout << PHWHERE << "caloTreeGen::process_event: Cluster " << m_clusterNode << " node is missing. Output related to this node will be empty" << std::endl;
    return 0;
  }

  // tower information
  TowerInfoContainer *emcTowerContainer;
  emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, m_emcTowerNode.c_str());
  if (!emcTowerContainer && storeEMCal)
  {
    std::cout << PHWHERE << "caloTreeGen::process_event: EMCal " << m_emcTowerNode << " node is missing. Output related to this node will be empty" << std::endl;
  }

  // grab all the towers and fill their energies.
  unsigned int tower_range = 0;
  if (storeEMCal && emcTowerContainer)
  {
    tower_range = emcTowerContainer->size();
    totalCaloEEMCal = 0;
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      unsigned int towerkey = emcTowerContainer->encode_key(iter);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
      double energy = emcTowerContainer->get_tower_at_channel(iter)->get_energy();
      if (energy < emcalThresh)
      {
        continue;
      }
      int time = emcTowerContainer->get_tower_at_channel(iter)->get_time();
      float chi2 = emcTowerContainer->get_tower_at_channel(iter)->get_chi2();
      float pedestal = emcTowerContainer->get_tower_at_channel(iter)->get_pedestal();
      totalCaloEEMCal += energy;
      m_emciPhi.push_back(iphi);
      m_emciEta.push_back(ieta);
      m_emcTowE.push_back(energy);
      m_emcTime.push_back(time);
      m_emcChi2.push_back(chi2);
      m_emcPed.push_back(pedestal);
    }
  }

  if (storeClusters && storeEMCal)
  {
    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter;
    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
    {
      RawCluster *recoCluster = clusterIter->second;

      CLHEP::Hep3Vector vertex(0, 0, 0);
      if (m_vertex != -9999)
      {
        vertex.setZ(m_vertex);
      }
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
      CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*recoCluster, vertex);

      float clusE = E_vec_cluster_Full.mag();
      if (clusE < clusterThresh)
      {
        continue;
      }

      float clusEcore = E_vec_cluster.mag();
      float clus_eta = E_vec_cluster.pseudoRapidity();
      float clus_phi = E_vec_cluster.phi();
      float clus_pt = E_vec_cluster.perp();
      float clus_chi = recoCluster->get_chi2();
      float nTowers = recoCluster->getNTowers();
      float maxTowerEnergy = getMaxTowerE(recoCluster);

      m_clusterE.push_back(clusE);
      m_clusterECore.push_back(clusEcore);
      m_clusterPhi.push_back(clus_phi);
      m_clusterEta.push_back(clus_eta);
      m_clusterPt.push_back(clus_pt);
      m_clusterChi.push_back(clus_chi);
      m_clusterNtow.push_back(nTowers);
      m_clusterTowMaxE.push_back(maxTowerEnergy);
      if (storeClusterDetails)
      {
        m_clusTowPhi.push_back(returnClusterTowPhi(recoCluster));
        m_clusTowEta.push_back(returnClusterTowEta(recoCluster));
        m_clusTowE.push_back(returnClusterTowE(recoCluster));
      }
    }
  }

  // tower information
  TowerInfoContainer *ohcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, m_ohcTowerNode);
  TowerInfoContainer *ihcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, m_ihcTowerNode.c_str());

  if ((!ohcTowerContainer || !ihcTowerContainer) && storeHCals)
  {
    std::cout << PHWHERE << "caloTreeGen::process_event: OHCal " << m_ohcTowerNode << " or IHCal " << m_ihcTowerNode << " node is missing. Output related to this node will be empty" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (storeHCals && ohcTowerContainer && ihcTowerContainer)
  {
    tower_range = ohcTowerContainer->size();
    totalCaloEOHCal = 0;
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      unsigned int towerkey = ohcTowerContainer->encode_key(iter);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
      int time = ohcTowerContainer->get_tower_at_channel(iter)->get_time();
      float chi2 = ohcTowerContainer->get_tower_at_channel(iter)->get_chi2();
      double energy = ohcTowerContainer->get_tower_at_channel(iter)->get_energy();
      if (energy < ohcalThresh)
      {
        continue;
      }
      float pedestal = ohcTowerContainer->get_tower_at_channel(iter)->get_pedestal();

      totalCaloEOHCal += energy;
      m_ohcTowE.push_back(energy);
      m_ohciTowEta.push_back(ieta);
      m_ohciTowPhi.push_back(iphi);
      m_ohcTime.push_back(time);
      m_ohcChi2.push_back(chi2);
      m_ohcPed.push_back(pedestal);
    }

    tower_range = ihcTowerContainer->size();
    totalCaloEIHCal = 0;
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      unsigned int towerkey = ihcTowerContainer->encode_key(iter);
      unsigned int ieta = TowerInfoDefs::getCaloTowerEtaBin(towerkey);
      unsigned int iphi = TowerInfoDefs::getCaloTowerPhiBin(towerkey);
      int time = ihcTowerContainer->get_tower_at_channel(iter)->get_time();
      float chi2 = ihcTowerContainer->get_tower_at_channel(iter)->get_chi2();
      double energy = ihcTowerContainer->get_tower_at_channel(iter)->get_energy();
      if (energy < ihcalThresh)
      {
        continue;
      }
      float pedestal = ihcTowerContainer->get_tower_at_channel(iter)->get_pedestal();

      totalCaloEIHCal += energy;
      m_ihcTowE.push_back(energy);
      m_ihciTowEta.push_back(ieta);
      m_ihciTowPhi.push_back(iphi);
      m_ihcTime.push_back(time);
      m_ihcChi2.push_back(chi2);
      m_ihcPed.push_back(pedestal);
    }
  }

  TowerInfoContainer *zdcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, m_zdcTowerNode.c_str());
  if (!zdcTowerContainer && storeZDC)
  {
    std::cout << PHWHERE << "caloTreeGen::process_event: ZDC " << m_zdcTowerNode << " node is missing. Output related to this node will be empty" << std::endl;
  }

  if (storeZDC && zdcTowerContainer)
  {
    tower_range = zdcTowerContainer->size();
    totalCaloEZDC = 0;
    for (unsigned int iter = 0; iter < tower_range; iter++)
    {
      if (iter < 16)
      {
        float energy = zdcTowerContainer->get_tower_at_channel(iter)->get_energy();
        m_zdcTowE.push_back(energy);
        unsigned int towerkey = zdcTowerContainer->encode_key(iter);
        unsigned int side = TowerInfoDefs::get_zdc_side(towerkey);
        m_zdcSide.push_back(side);
        totalCaloEZDC += energy;
      }
      if (iter > 15 && iter < 48)
      {
        // smd north stuff
        float energy = zdcTowerContainer->get_tower_at_channel(iter)->get_energy();
        m_smdE.push_back(energy);
        unsigned int towerkey = zdcTowerContainer->encode_key(iter);
        unsigned int side = TowerInfoDefs::get_zdc_side(towerkey);
        m_smdSide.push_back(side);
      }
    }
  }

  Gl1Packet *gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, m_trigNode.c_str());
  if (!gl1PacketInfo && storeTrig)
  {
    std::cout << PHWHERE << "caloTreeGen::process_event: Gl1 " << m_trigNode << " node is missing. Output related to this node will be empty" << std::endl;
  }

  if (storeTrig && gl1PacketInfo)
  {
    uint64_t triggervec = gl1PacketInfo->getScaledVector();
    for (int i = 0; i < 64; i++)
    {
      bool trig_decision = ((triggervec & 0x1U) == 0x1U);
      m_triggerVector.push_back(trig_decision);
      triggervec = (triggervec >> 1U) & 0xffffffffU;
    }
  }
  T->Fill();
  eventNum++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::ResetEvent(PHCompositeNode * /*topNode*/)
{
  m_clusterE.clear();
  m_clusterPhi.clear();
  m_clusterEta.clear();
  m_clusterPt.clear();
  m_clusterChi.clear();
  m_clusterTowMaxE.clear();
  m_clusterNtow.clear();
  m_clusterECore.clear();

  m_emcTowE.clear();
  m_emciEta.clear();
  m_emciPhi.clear();
  m_emcTime.clear();
  m_emcChi2.clear();
  m_emcPed.clear();

  m_ihcTowE.clear();
  m_ihciTowEta.clear();
  m_ihciTowPhi.clear();
  m_ihcTime.clear();
  m_ihcChi2.clear();
  m_ihcPed.clear();

  m_ohcTowE.clear();
  m_ohciTowEta.clear();
  m_ohciTowPhi.clear();
  m_ohcTime.clear();
  m_ohcChi2.clear();
  m_ohcPed.clear();

  m_clusTowPhi.clear();
  m_clusTowEta.clear();
  m_clusTowE.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int caloTreeGen::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "caloTreeGen::End(PHCompositeNode *topNode) Saving TTree" << std::endl;

  out->cd();
  T->Write();
  zVertex->Write();
  out->Close();
  delete out;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
float caloTreeGen::getMaxTowerE(RawCluster *cluster)
{
  RawCluster::TowerConstRange towers = cluster->get_towers();
  RawCluster::TowerConstIterator toweriter;

  float maxEnergy = 0;
  for (toweriter = towers.first; toweriter != towers.second; toweriter++)
  {
    float towE = toweriter->second;

    if (towE > maxEnergy)
    {
      maxEnergy = towE;
    }
  }
  return maxEnergy;
}
//____________________________________________________________________________..
std::vector<int> caloTreeGen::returnClusterTowEta(RawCluster *cluster)
{
  RawCluster::TowerConstRange towers = cluster->get_towers();
  RawCluster::TowerConstIterator toweriter;

  std::vector<int> towerIDsEta;
  for (toweriter = towers.first; toweriter != towers.second; toweriter++)
  {
    towerIDsEta.push_back(RawTowerDefs::decode_index1(toweriter->first));
  }

  return towerIDsEta;
}
//____________________________________________________________________________..
std::vector<int> caloTreeGen::returnClusterTowPhi(RawCluster *cluster)
{
  RawCluster::TowerConstRange towers = cluster->get_towers();
  RawCluster::TowerConstIterator toweriter;

  std::vector<int> towerIDsPhi;
  for (toweriter = towers.first; toweriter != towers.second; toweriter++)
  {
    towerIDsPhi.push_back(RawTowerDefs::decode_index2(toweriter->first));
  }
  return towerIDsPhi;
}
//____________________________________________________________________________..
std::vector<float> caloTreeGen::returnClusterTowE(RawCluster *cluster)
{
  RawCluster::TowerConstRange towers = cluster->get_towers();
  RawCluster::TowerConstIterator toweriter;

  std::vector<float> towerE;
  for (toweriter = towers.first; toweriter != towers.second; toweriter++)
  {
    towerE.push_back(toweriter->second);
  }

  return towerE;
}
