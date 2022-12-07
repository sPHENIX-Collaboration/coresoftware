#include "CaloCalibEmc_Pi0.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector

#include <algorithm>  // for max, max_element
#include <cmath>      // for abs
#include <cstdlib>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>   // for vector

//using namespace std;

//____________________________________________________________________________..
CaloCalibEmc_Pi0::CaloCalibEmc_Pi0(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_Filename(filename)
{
  for (int nj = 0; nj < 96; nj++)
  {
    eta_hist[nj] = nullptr;
    for (int nk = 0; nk < 256; nk++)
    {
      cemc_hist_eta_phi[nj][nk] = nullptr;
    }
  }
}

//____________________________________________________________________________..
int CaloCalibEmc_Pi0::InitRun(PHCompositeNode *topNode)
{
  std::cout << "LiteCaloEval::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  m_ievent = 0;

  cal_output = new TFile(m_Filename.c_str(), "RECREATE");

  pairInvMassTotal = new TH1F("pairInvMassTotal", "Pair Mass Histo", 70, 0.0, 0.7);
  mass_eta = new TH2F("mass_eta", "2d Pair Mass Histo", 70, 0.0, 0.7, 400, -1.5, 1.5);
  mass_eta_phi = new TH3F("mass_eta_phi", "3d Pair Mass Histo", 70, 0.0, 0.7, 150, -1.5, 1.5, 256, -3.142, 3.142);

  // histo to record every tower by tower locations
  for (int i = 0; i < 96; i++)  // eta rows
  {
    for (int j = 0; j < 258; j++)  // phi columns
    {
      TString i1;
      TString j1;
      i1.Form("%d", i);
      j1.Form("%d", j);
      TString hist_name = "emc_ieta" + i1 + "_phi" + j1;

      cemc_hist_eta_phi[i][j] = new TH1F(hist_name.Data(), "Hist_ieta_phi_", 70, 0.0, 0.7);
    }
  }

  // histo to record each eta locations (with all phis included in each)
  for (int i = 0; i < 96; i++)
  {
    gStyle->SetOptFit(1);
    TString a;
    a.Form("%d", i);
    TString b = "eta_" + a;

    eta_hist[i] = new TH1F(b.Data(), "eta and all phi's", 70, 0.0, 0.7);
  }

  if (topNode != 0)
  {
    // TTree declare
    _eventTree = new TTree("_eventTree", "An event level info Tree");
    // TTree branches
    _eventTree->Branch("_eventNumber", &_eventNumber, "_eventNumber/I");
    _eventTree->Branch("_nClusters", &_nClusters, "_nClusters/I");
    _eventTree->Branch("_clusterIDs", _clusterIDs, "_clusterIDs[_nClusters]/I");
    _eventTree->Branch("_clusterEnergies", _clusterEnergies, "_clusterEnergies[_nClusters]/F");
    _eventTree->Branch("_clusterPts", _clusterPts, "_clusterPts[_nClusters]/F");
    _eventTree->Branch("_clusterEtas", _clusterEtas, "_clusterEtas[_nClusters]/I");
    _eventTree->Branch("_clusterPhis", _clusterPhis, "_clusterPhis[_nClusters]/I");
    _eventTree->Branch("_maxTowerEtas", _maxTowerEtas, "_maxTowerEtas[_nClusters]/I");
    _eventTree->Branch("_maxTowerPhis", _maxTowerPhis, "_maxTowerPhis[_nClusters]/I");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloCalibEmc_Pi0::process_event(PHCompositeNode *topNode)
{
  if (m_ievent % 50 == 0)
  {
    std::cout << std::endl;
    std::cout << "Beginning of the event " << m_ievent << std::endl;
    std::cout << "====================================" << std::endl;
  }

  _eventNumber = m_ievent++;

  // create a cluster object
  RawClusterContainer *recal_clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_POS_COR_CEMC");

  // create a tower object
  std::string towernode = "TOWER_CALIB_" + _caloname;
  RawTowerContainer *_towers = findNode::getClass<RawTowerContainer>(topNode, towernode.c_str());
  if (!_towers)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
  }

  // create a tower geometry object
  std::string towergeomnode = "TOWERGEOM_" + _caloname;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode.c_str());
  if (!towergeom)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnode << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get Vertex
  float vx = 0;
  float vy = 0;
  float vz = 0;
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (vertexmap)
  {
    if (!vertexmap->empty())
    {
      GlobalVertex *vtx = (vertexmap->begin()->second);
      vx = vtx->get_x();
      vy = vtx->get_y();
      vz = vtx->get_z();
    }
  }
  /////////////////////////////

  //	std::cout << "vtx " << vx  << "  " << vy << std::endl;
  // CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
  CLHEP::Hep3Vector vertex(vx, vy, vz);

  // ------------------------------

  // loop over the clusters
  RawClusterContainer::ConstRange t_rbegin_end = recal_clusters->getClusters();
  RawClusterContainer::ConstIterator t_rclusiter;

  RawCluster *savCs[10000];  // savingClusters that has 1 GeV or more
  int iCs = 0;
  int inCs = 0;

  // saving the clusters
  for (t_rclusiter = t_rbegin_end.first; t_rclusiter != t_rbegin_end.second; ++t_rclusiter)
  {
    RawCluster *t_recalcluster = t_rclusiter->second;

    float cluse = t_recalcluster->get_energy();
    if (cluse > 0.1) inCs++;

    if (cluse > 1.0) savCs[iCs++] = t_recalcluster;
  }

  _nClusters = iCs;
  if (_nClusters > 60)
    return Fun4AllReturnCodes::EVENT_OK;

  // looping on the saved clusters savCs[]
  // outer loop (we want to do pair of the loops)
  for (int jCs = 0; jCs < iCs; jCs++)
  {
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetEVec(*savCs[jCs], vertex);

    // CLHEP::Hep3Vector E_vec_cluster(PxCluster, PyCluster, ...);

    _clusterIDs[jCs] = savCs[jCs]->get_id();
    _clusterEnergies[jCs] = savCs[jCs]->get_energy();

    float tt_clus_energy = E_vec_cluster.mag();
    if (tt_clus_energy > 1.0)  // it is greater than 1.0 no need to do if here but lets keep tho
    {
      float tt_clus_eta = E_vec_cluster.pseudoRapidity();
      float tt_clus_pt = E_vec_cluster.perp();
      float tt_clus_phi = E_vec_cluster.getPhi();
      _clusterPts[jCs] = tt_clus_pt;
      _clusterEtas[jCs] = tt_clus_eta;
      _clusterPhis[jCs] = tt_clus_phi;

      // another loop to go into the saved cluster
      for (int kCs = 0; kCs < iCs; kCs++)
      {
        if (jCs == kCs) continue;

        CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetEVec(*savCs[kCs], vertex);

        float tt2_clus_energy = E_vec_cluster2.mag();
        if (tt2_clus_energy > 1.0)  // again this will be greater than 1.0, but lets keep
        {
          // lets do alpha cut here: this is needed tho
          alphaCut = abs(tt_clus_energy - tt2_clus_energy) / (tt_clus_energy + tt_clus_energy);
          if (alphaCut <= 0.8)
          {
            float tt2_clus_eta = E_vec_cluster2.pseudoRapidity();
            float tt2_clus_pt = E_vec_cluster2.perp();
            float tt2_clus_phi = E_vec_cluster2.getPhi();

            TLorentzVector pho1, pho2, pi0lv;

            pho1.SetPtEtaPhiE(tt_clus_pt, tt_clus_eta, tt_clus_phi, tt_clus_energy);
            pho2.SetPtEtaPhiE(tt2_clus_pt, tt2_clus_eta, tt2_clus_phi, tt2_clus_energy);

            if (pho1.DeltaR(pho2) > 0.8) continue;

            pi0lv = pho1 + pho2;
            float pairInvMass = pi0lv.M();

            //vector to hold all the towers etas, phis, and energy in this cluster
            std::vector<int> toweretas;
            std::vector<int> towerphis;
            std::vector<float> towerenergies;

            // loop over the towers from the outer loop cluster
            // and find the max tower location and save the
            // histogram on that max tower location with this
            // invariant mass

            RawCluster::TowerConstRange towers = savCs[jCs]->get_towers();
            RawCluster::TowerConstIterator toweriter;

            for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
            {
              RawTower *tower = _towers->getTower(toweriter->first);

              int towereta = tower->get_bineta();
              int towerphi = tower->get_binphi();
              double towerenergy = tower->get_energy();

              // put the eta, phi, energy into corresponding vectors
              toweretas.push_back(towereta);
              towerphis.push_back(towerphi);
              towerenergies.push_back(towerenergy);
            }

            // std::cout << std::endl;
            // std::cout << "Cluster energy: " << tt_clus_energy << std::endl;
            // std::cout << "Total number of towers (getNTowers()): " << savCs[jCs]->getNTowers() << std::endl;
            // std::cout << "Total number of towers size(toweretas): " << toweretas.size() << std::endl;
            // float maxTowerEnergy = *max_element(towerenergies.begin(), towerenergies.end());
            // std::cout << "The maxTowerEnergy: " << maxTowerEnergy << std::endl;
            int maxTowerIndex = max_element(towerenergies.begin(), towerenergies.end()) - towerenergies.begin();
            maxTowerEta = toweretas[maxTowerIndex];
            maxTowerPhi = towerphis[maxTowerIndex];

            _maxTowerEtas[jCs] = maxTowerEta;
            _maxTowerPhis[jCs] = maxTowerPhi;

            if (tt_clus_energy > 2.5 && tt2_clus_energy > 1.5)
            {
              pairInvMassTotal->Fill(pairInvMass);
              mass_eta->Fill(pairInvMass, tt_clus_eta);
              mass_eta_phi->Fill(pairInvMass, tt_clus_eta, tt_clus_phi);
            }

            // fill the tower by tower histograms with invariant mass
            cemc_hist_eta_phi[maxTowerEta][maxTowerPhi]->Fill(pairInvMass);
            eta_hist[maxTowerEta]->Fill(pairInvMass);
          }
        }
      }
    }
  }
  _eventTree->Fill();

  //    m_ievent++;

  // pause few seconds
  //	std::chrono::seconds dura(2);
  //	std::this_thread::sleep_for(dura);
  //	std::std::cout << "Waited 2 sec\n";

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloCalibEmc_Pi0::End(PHCompositeNode * /*topNode*/)
{
  cal_output->cd();
  //	_eventTree->Write();
  cal_output->Write();
  cal_output->Close();
  delete cal_output;

  return Fun4AllReturnCodes::EVENT_OK;
}

//______________________________________________________________________________..
void CaloCalibEmc_Pi0::Loop(TString filename, int nevts)
{
  TFile *f = new TFile(filename);

  TTree *t1 = (TTree *) f->Get("_eventTree");

  // Declaration of leaf types
  Int_t eventNumber;
  Int_t nClusters;
  Int_t clusterIDs[10000];         //[_nClusters]
  Float_t clusterEnergies[10000];  //[_nClusters]
  Float_t clusterPts[10000];
  Int_t clusterEtas[10000];  //[_nClusters]
  Int_t clusterPhis[10000];  //[_nClusters]
  Int_t maxTowerEtas[10000];
  Int_t maxTowerPhis[10000];

  // Set Branches
  t1->SetBranchAddress("_eventNumber", &eventNumber);
  t1->SetBranchAddress("_nClusters", &nClusters);
  t1->SetBranchAddress("_clusterIDs", clusterIDs);
  t1->SetBranchAddress("_clusterEnergies", clusterEnergies);
  t1->SetBranchAddress("_clusterPts", clusterPts);
  t1->SetBranchAddress("_clusterEtas", clusterEtas);
  t1->SetBranchAddress("_clusterPhis", clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", maxTowerPhis);

  // pre-loop to save all the clusters LorentzVector

  TLorentzVector *savClusLV[10000];

  int nEntries = (int) t1->GetEntriesFast();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
    nevts2 = nEntries;

  for (int i = 0; i < nevts2; i++)
  {
    // load the ith instance of the TTree
    t1->GetEntry(i);

    // calibration correction will be applied here
    //
    int nClusters_tmp = _nClusters;
    for (int j = 0; j < nClusters_tmp; j++)
    {
      // float px, py, pz;
      float pt, eta, phi, E;
      pt = _clusterPts[j];
      eta = _clusterEtas[j];
      phi = _clusterPhis[j];
      E = _clusterEnergies[j];
      // px  =  pt * cos(phi);
      // py  =  pt * sin(phi);
      // pz  =  pt * sinh(eta);
      savClusLV[j] = new TLorentzVector();
      savClusLV[j]->SetPtEtaPhiE(pt, eta, phi, E);
    }

    // Next we need to have the clusterCut to only take
    // certain clusters of certain energy
    // But lets keep it simple for now

    TLorentzVector *pho1, *pho2;
    int iCs = nClusters_tmp;
    for (int jCs = 0; jCs < iCs; jCs++)
    {
      pho1 = savClusLV[jCs];

      // another loop to go into the saved cluster
      for (int kCs = 0; kCs < iCs; kCs++)
      {
        if (jCs == kCs) continue;

        pho2 = savClusLV[kCs];

        TLorentzVector pi0lv;

        if (pho1->DeltaR(*pho2) > 0.8) continue;

        pi0lv = *pho1 + *pho2;
        float pairInvMass = pi0lv.M();
        //std::cout << pairInvMass << std::endl;

        // fill the tower by tower histograms with invariant mass
        cemc_hist_eta_phi[_maxTowerEtas[jCs]][_maxTowerPhis[jCs]]->Fill(pairInvMass);
        eta_hist[_maxTowerEtas[jCs]]->Fill(pairInvMass);
      }
    }
  }
}

// _______________________________________________________________..
void CaloCalibEmc_Pi0::FittingHistos()
{
  TFile *parFile = new TFile("parFile", "RECREATE");

  TTree *parTree = new TTree("parTree", "Tree with fit mean saved");

  int ieta, iphi;
  float corrval[100000];

  parTree->Branch("ieta", &ieta, "ieta/I");
  parTree->Branch("iphi", &iphi, "iphi/I");
  parTree->Branch("corrval", corrval, "corrval[ieta*1000+iphi]");

  TF1 *fit_result;

  for (int i = 0; i < 96; i++)
  {
    for (int j = 0; j <= 258; j++)
    {
      if (j != 258)
      {
        cemc_hist_eta_phi[i][j]->Fit("gaus", "", "", 0.095, 0.175);
        fit_result = cemc_hist_eta_phi[i][j]->GetFunction("gaus");
        ieta = i;
        iphi = j;
        corrval[i * 1000 + j] = fit_result->GetParameter(1);
      }
      else
      {
        eta_hist[i]->Fit("gaus", "", "", 0.095, 0.175);
        fit_result = eta_hist[i]->GetFunction("gaus");
        ieta = i;
        iphi = j;
        corrval[i * 1000 + j] = fit_result->GetParameter(1);
      }
    }
  }
  parTree->Fill();
  parTree->Write();
  parFile->Close();
}
