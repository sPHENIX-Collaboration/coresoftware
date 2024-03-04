#include "CaloCalibEmc_Pi0.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector

#include <algorithm>  // for max, max_element
#include <cmath>      // for abs
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <string>   // for string
#include <utility>  // for pair
#include <vector>   // for vector

//____________________________________________________________________________..
CaloCalibEmc_Pi0::CaloCalibEmc_Pi0(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_Filename(filename)
  , m_cent_nclus_cut(350)
{
  eta_hist.fill(nullptr);
  // the 2d filling took a bit of googling but it works:
  std::for_each(cemc_hist_eta_phi.begin(), cemc_hist_eta_phi.end(), [](auto &row)
                { row.fill(nullptr); });
}

//____________________________________________________________________________..
int CaloCalibEmc_Pi0::InitRun(PHCompositeNode *topNode)
{
  std::cout << "LiteCaloEval::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  m_ievent = 0;

  cal_output = new TFile(m_Filename.c_str(), "RECREATE");

  pairInvMassTotal = new TH1F("pairInvMassTotal", "Pair Mass Histo", 130, 0.0, 1.3);
  mass_eta = new TH2F("mass_eta", "2d Pair Mass Histo", 70, 0.0, 0.7, 400, -1.5, 1.5);
  mass_eta_phi = new TH3F("mass_eta_phi", "3d Pair Mass Histo", 70, 0.0, 0.7, 150, -1.5, 1.5, 256, -3.142, 3.142);
  h_totalClusters = new TH1F("h_totalClusters", "Histogram of total number of clusters", 600, 0, 600);
  pt1_ptpi0_alpha = new TH3F("pt1_ptpi0_alpha", "Photon1 pT, pi0 pT, alpha", 40, 0., 4., 40, 0., 4., 10, 0., 1.);
  fitp1_eta_phi2d = new TH2F("fitp1_eta_phi2d", "fit p1 eta phi", 16, 0, 16, 16, 0, 16);

  // temporarily don't need these till Run2 and they take up a lot of space
  //  // histo to record every tower by tower locations
  //  for (int i = 0; i < 96; i++)  // eta rows
  //  {
  //    for (int j = 0; j < 258; j++)  // phi columns
  //    {
  //      std::string hist_name = std::string("emc_ieta") + std::tostring(i)
  //          + std::string("_phi") + std::to_string(j);

  //     cemc_hist_eta_phi[i][j] = new TH1F(hist_name.c_str()(), "Hist_ieta_phi_", 70, 0.0, 0.7);
  //   }
  // }

  // histo to record each eta locations (with all phis included in each)
  for (int i = 0; i < 96; i++)
  {
    gStyle->SetOptFit(1);
    std::string b = std::string("eta_") + std::to_string(i);

    eta_hist.at(i) = new TH1F(b.c_str(), "eta and all phi's", 130, 0.0, 1.3);
  }

  if (topNode != nullptr)
  {
    // TTree declare
    _eventTree = new TTree("_eventTree", "An event level info Tree");
    // TTree branches
    _eventTree->Branch("_eventNumber", &_eventNumber, "_eventNumber/I");
    _eventTree->Branch("_nClusters", &_nClusters, "_nClusters/I");
    //_eventTree->Branch("_nTowers", &_nTowers, "_nTowers[_nClusters]/I");
    _eventTree->Branch("_clusterIDs", _clusterIDs, "_clusterIDs[_nClusters]/I");
    _eventTree->Branch("_clusterEnergies", _clusterEnergies, "_clusterEnergies[_nClusters]/F");
    _eventTree->Branch("_clusterPts", _clusterPts, "_clusterPts[_nClusters]/F");
    _eventTree->Branch("_clusterEtas", _clusterEtas, "_clusterEtas[_nClusters]/F");
    _eventTree->Branch("_clusterPhis", _clusterPhis, "_clusterPhis[_nClusters]/F");
    _eventTree->Branch("_maxTowerEtas", _maxTowerEtas, "_maxTowerEtas[_nClusters]/I");
    _eventTree->Branch("_maxTowerPhis", _maxTowerPhis, "_maxTowerPhis[_nClusters]/I");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloCalibEmc_Pi0::process_event(PHCompositeNode *topNode)
{
  if (m_ievent % 50 == 0)
  //  if (true)
  {
    std::cout << std::endl;
    std::cout << "Beginning of the event " << m_ievent << std::endl;
    std::cout << "====================================" << std::endl;
  }

  _eventNumber = m_ievent++;

  std::string clusnodename = "CLUSTER_CEMC";
  if (m_UseTowerInfo && !_inputnodename.empty())
  {
    //    clusnodename = "CLUSTERINFO_CEMC";
    //    clusnodename = "CLUSTERINFO2_CEMC";
    clusnodename = _inputnodename;
  }

  // create a cluster object
  RawClusterContainer *recal_clusters = findNode::getClass<RawClusterContainer>(topNode, clusnodename);

  if (!recal_clusters)
  {
    std::cout << PHWHERE << "EMCal cluster node is missing, can't collect EMCal clusters" << std::endl;
  }

  // create a tower object
  std::string towernode = "TOWER_CALIB_" + _caloname;
  RawTowerContainer *_towers = findNode::getClass<RawTowerContainer>(topNode, towernode);
  if (!_towers && m_UseTowerInfo < 1)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
  }

  // create a tower geometry object
  std::string towergeomnode = "TOWERGEOM_" + _caloname;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode);
  if (!towergeom && m_UseTowerInfo < 1)
  {
    std::cout << PHWHERE << ": Could not find node " << towergeomnode << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // if using towerinfo
  //  create a tower object
  //   towernode = "TOWERINFO_CALIB_" + _caloname;
  towernode = "TOWERS_Calib_" + _caloname;
  if (!_inputtownodename.empty())
  {
    towernode = _inputtownodename;
  }

  TowerInfoContainer *_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towernode);
  if (!_towerinfos && m_UseTowerInfo > 0)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
  }

  if (m_ievent < 3)
  {
    std::cout << "using these nodes for cluster/tower: " << clusnodename << " " << towernode << std::endl;
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

    float cluse = t_recalcluster->get_ecore();
    //    std::cout << "clus " << cluse << std::endl;

    if (cluse > 0.5)
    {
      inCs++;
    }

    if (cluse > 0.6 && t_recalcluster->get_chi2() < 4)
    {
      savCs[iCs++] = t_recalcluster;
    }
  }

  _nClusters = iCs;

  if (_nClusters > m_cent_nclus_cut)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // looping on the saved clusters savCs[]
  // outer loop (we want to do pair of the loops)
  for (int jCs = 0; jCs < iCs; jCs++)
  {
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*savCs[jCs], vertex);

    // CLHEP::Hep3Vector E_vec_cluster(PxCluster, PyCluster, ...);

    _clusterIDs[jCs] = savCs[jCs]->get_id();

    // vector to hold all the towers etas, phis, and energy in this cluster
    std::vector<int> toweretas;
    std::vector<int> towerphis;
    std::vector<float> towerenergies;

    // loop over the towers from the outer loop cluster
    // and find the max tower location and save the
    // histogram on that max tower location with this
    // invariant mass

    RawCluster::TowerConstRange towers = savCs[jCs]->get_towers();
    RawCluster::TowerConstIterator toweriter;

    int towerieta = -1;
    int toweriphi = -1;

    for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
    {
      if (m_UseTowerInfo < 1)
      {
        RawTower *tower = _towers->getTower(toweriter->first);
        towerieta = tower->get_bineta();
        toweriphi = tower->get_binphi();
      }
      else
      {
        // this probably works for case above as well but to be checked...
        // (ie probably  no need for RawTower pointer or the if statement)
        toweriphi = RawTowerDefs::decode_index2(toweriter->first);  // index2 is phi in CYL
        towerieta = RawTowerDefs::decode_index1(toweriter->first);  // index1 is eta in CYL
      }

      double towerenergy = toweriter->second;

      // put the eta, phi, energy into corresponding vectors
      toweretas.push_back(towerieta);
      towerphis.push_back(toweriphi);
      towerenergies.push_back(towerenergy);
    }

    // cout << endl;
    // cout << "Cluster energy: " << tt_clus_energy << endl;
    // cout << "Total number of towers (getNTowers()): " << savCs[jCs]->getNTowers() << endl;
    // cout << "Total number of towers size(toweretas): " << toweretas.size() << endl;
    // float maxTowerEnergy = *max_element(towerenergies.begin(), towerenergies.end());
    // cout << "The maxTowerEnergy: " << maxTowerEnergy << endl;

    int maxTowerIndex = max_element(towerenergies.begin(), towerenergies.end()) - towerenergies.begin();
    maxTowerEta = toweretas[maxTowerIndex];
    maxTowerPhi = towerphis[maxTowerIndex];

    _maxTowerEtas[jCs] = maxTowerEta;
    _maxTowerPhis[jCs] = maxTowerPhi;

    float tt_clus_energy = E_vec_cluster.mag();
    float tt_clus_eta = E_vec_cluster.pseudoRapidity();
    float tt_clus_pt = E_vec_cluster.perp();
    float tt_clus_phi = E_vec_cluster.getPhi();

    // if (tt_clus_energy > 0.029)
    // 	{

    // 	if (m_UseTowerInfo < 1)
    // 	  std::cout <<  "rawtJF clus " << jCs << " " << tt_clus_energy <<  " phi: "
    // 		    << tt_clus_phi << " eta: " << tt_clus_eta <<  " ieta  "
    // 		    << maxTowerEta << " iphi " << maxTowerPhi
    // 		    << " towsize " << towerenergies.size() << std::endl;
    // 	else
    // 	  std::cout <<  "infoJF clus " << jCs << " " << tt_clus_energy <<  " phi: "
    // 		    << tt_clus_phi << " eta: " << tt_clus_eta <<  " ieta  "
    // 		    << maxTowerEta << " iphi " << maxTowerPhi
    // 		    << " towsize " << towerenergies.size() << std::endl;

    // 	}

    _clusterEnergies[jCs] = tt_clus_energy;
    _clusterPts[jCs] = tt_clus_pt;
    _clusterEtas[jCs] = tt_clus_eta;
    _clusterPhis[jCs] = tt_clus_phi;

    if (tt_clus_pt > 0.9 && iCs > 400000)  // while we are using the ntuples we don't need the
    // rest of this fn --but this will likely change soon.
    {
      // another loop to go into the saved cluster
      for (int kCs = 0; kCs < iCs; kCs++)
      {
        if (jCs == kCs)
        {
          continue;
        }

        CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*savCs[kCs], vertex);

        float tt2_clus_energy = E_vec_cluster2.mag();
        if (tt2_clus_energy > 0.6)  // again this will be greater than 1.0, but lets keep
        {
          // lets do alpha cut here: this is needed tho
          alphaCut = std::abs(tt_clus_energy - tt2_clus_energy) / (tt_clus_energy + tt_clus_energy);
          if (alphaCut <= 0.5)
          {
            float tt2_clus_eta = E_vec_cluster2.pseudoRapidity();
            float tt2_clus_pt = E_vec_cluster2.perp();
            float tt2_clus_phi = E_vec_cluster2.getPhi();

            TLorentzVector pho1, pho2, pi0lv;

            pho1.SetPtEtaPhiE(tt_clus_pt, tt_clus_eta, tt_clus_phi, tt_clus_energy);
            pho2.SetPtEtaPhiE(tt2_clus_pt, tt2_clus_eta, tt2_clus_phi, tt2_clus_energy);

            if (pho1.DeltaR(pho2) > 0.8)
            {
              continue;
            }
            // if (pho1.Eta()/pho2.Eta() < 0) continue;

            pi0lv = pho1 + pho2;
            float pairInvMass = pi0lv.M();

            pt1_ptpi0_alpha->Fill(tt_clus_pt, fabs(pi0lv.Pt()), alphaCut);

            //		      if (tt_clus_energy > 1.3 && tt2_clus_energy > 1.0 && fabs(pi0lv.Pt()) > 1.8)
            if (fabs(pi0lv.Pt()) > 1.0)
            {
              pairInvMassTotal->Fill(pairInvMass);
              mass_eta->Fill(pairInvMass, tt_clus_eta);
              //  mass_eta_phi->Fill(pairInvMass, tt_clus_eta, tt_clus_phi);

              // fill the tower by tower histograms with invariant mass
              cemc_hist_eta_phi.at(maxTowerEta).at(maxTowerPhi)->Fill(pairInvMass);
              eta_hist.at(maxTowerEta)->Fill(pairInvMass);
            }
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
int CaloCalibEmc_Pi0::End(PHCompositeNode *topNode)
{
  if (topNode == nullptr && f_temp)
  {
    cal_output->Close();
    f_temp->Close();
    delete f_temp;
    delete cal_output;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  cal_output->cd();
  //	_eventTree->Write();
  cal_output->Write();
  cal_output->Close();
  delete cal_output;

  return Fun4AllReturnCodes::EVENT_OK;
}

//______________________________________________________________________________..
void CaloCalibEmc_Pi0::Loop(int nevts, const std::string &filename, TTree *intree, const std::string &incorrFile)
{
  // KINEMATIC CUTS ON PI0's ARE LISTED BELOW IN ONE SECTION OF COMMENTS
  //
  //  search for  "CUTS FOLLOW HERE"
  //
  //  they are centrality dependent (centrality = nclusters > minpt)  so they can only
  //  defined after loading up this values from the ntuples
  //  so they can't be moved here
  //

  //  std::arrays have their indices backward, this is the old float myaggcorr[96][260];
  std::array<std::array<float, 260>, 96> myaggcorr{};
  std::for_each(myaggcorr.begin(), myaggcorr.end(), [](auto &row)
                { row.fill(1.); });

  std::cout << "running w/ corr file? : " << incorrFile << std::endl;

  if (!incorrFile.empty())
  {
    TFile *infileNt = new TFile(incorrFile.c_str());
    std::cout << "loaded incorrFile " << infileNt << std::endl;

    float myieta;
    float myiphi;
    float mycorr;
    float myaggcv;

    TNtuple *innt_corrVals = (TNtuple *) infileNt->Get("nt_corrVals");

    innt_corrVals->SetBranchAddress("tower_eta", &myieta);
    innt_corrVals->SetBranchAddress("tower_phi", &myiphi);
    innt_corrVals->SetBranchAddress("corr_val", &mycorr);
    innt_corrVals->SetBranchAddress("agg_cv", &myaggcv);

    int ntCorrs = innt_corrVals->GetEntries();

    for (int ij = 0; ij < ntCorrs; ij++)
    {
      innt_corrVals->GetEntry(ij);
      int ci = (int) myieta;
      int cj = (int) myiphi;
      myaggcorr.at(ci).at(cj) = myaggcv;
      if (ij > ntCorrs - 2 || ij == ntCorrs / 2)
      {
        std::cout << "loaded corrs eta,phi,aggcv " << myieta
                  << " " << myiphi << " " << myaggcv << std::endl;
      }
    }

    infileNt->Close();
    delete infileNt;
  }

  std::cout << "in loop" << std::endl;

  TTree *t1 = intree;
  if (!intree)
  {
    TFile *f = new TFile(filename.c_str());
    t1 = (TTree *) f->Get("_eventTree");
  }

  // Set Branches
  //  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  //  t1->SetBranchAddress("_clusterIDs", _clusterIDs);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);

  // pre-loop to save all the clusters LorentzVector

  TLorentzVector *savClusLV[10000];

  //  int nEntries = (int) t1->GetEntriesFast();
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
  {
    nevts2 = nEntries;
  }

  // keeping track of discarded clusters for v7
  int discarded_clusters = 0;

  for (int i = 0; i < nevts2; i++)
  {
    // load the ith instance of the TTree
    t1->GetEntry(i);

    if ((i % 10 == 0 && i < 200) || (i % 100 == 0 && i < 1000) || (i % 1000 == 0 && i < 37003) || i % 10000 == 0)
    {
      std::cout << "evt no " << i << std::endl;
    }
    // calibration correction will be applied here

    int nClusters = _nClusters;

    // see below this is like centrality cut, but currently need central events
    // as well as peripheral to maximize statistical power
    if (nClusters > 1000)
    {
      discarded_clusters += 1;
      continue;
    }

    float pt1cut = 0, pt2cut = 0;

    for (int j = 0; j < nClusters; j++)
    {
      // float px, py, pz;
      float pt, eta, phi, E, aggcv;
      pt = _clusterPts[j];
      eta = _clusterEtas[j];
      phi = _clusterPhis[j];
      E = _clusterEnergies[j];
      // px  =  pt * cos(phi);
      // py  =  pt * sin(phi);
      // pz  =  pt * sinh(eta);
      // pt *= myaggcorr[
      aggcv = myaggcorr.at(_maxTowerEtas[j]).at(_maxTowerPhis[j]);

      // std::cout << "aggcv applied: " << aggcv << std::endl;

      // eta slice shifts test
      //   int ket = _maxTowerEtas[j]/4;
      //   int jket = ket %4;
      //  if ((ket/4)%2==1)
      //	jket = 4-ket%4;
      //   int pjj = _maxTowerEtas[j]%4 - 1;
      //       aggcv *= 0.86+jket*0.11 + 0.02*pjj;

      pt *= aggcv;
      E *= aggcv;

      savClusLV[j] = new TLorentzVector();
      savClusLV[j]->SetPtEtaPhiE(pt, eta, phi, E);
    }

    TLorentzVector *pho1, *pho2;
    int iCs = nClusters;
    for (int jCs = 0; jCs < iCs; jCs++)
    {
      pho1 = savClusLV[jCs];
      /////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////
      // *********************************
      //
      //  CUTS FOLLOW HERE (e.g. pt cuts)
      //
      //*************************************
      ///////////////////////////////////

      // centrality dependent pt cuts designed to keep
      // statistical cluster count   contribution (& sig/bkg)
      // constant with all centrality
      // in order to maximize statistical power i.e. using all events
      // in the calibration not just peripheral events.
      // this is neccessary for the summer 23 data because
      // the event rate was small and the total statistics per
      // stable calibration period (typically a daq run-length) is small

      float modCutFactor = 1.0;

      if (iCs < 30)
      {
        // pt1cut =  1.65*modCutFactor;
        // pt2cut  = 0.8*modCutFactor;

        pt1cut = 1.3 * modCutFactor;
        pt2cut = 0.7 * modCutFactor;
      }
      else
      {
        // pt1cut = 1.65*modCutFactor +  1.4*(iCs-29)/200.0*modCutFactor;
        // pt2cut = 0.8*modCutFactor +  1.4*(iCs-29)/200.0*modCutFactor;

        pt1cut = 1.3 * modCutFactor + 1.4 * (iCs - 29) / 200.0 * modCutFactor;
        pt2cut = 0.7 * modCutFactor + 1.4 * (iCs - 29) / 200.0 * modCutFactor;
      }

      float pi0ptcut = 1.22 * (pt1cut + pt2cut);

      // energy asymmetry alpha cut
      float alphacutval = 0.6;

      float deltaRconecut = 1.1;  // 2-gamma opening angle(dR) cut
      // value relevant for background extent in mass
      //  not in peak area.

      ////////////////////////////////////////////////////////
      //////////////////////////////////
      //   END CUTS
      ///////////////////////////////////////
      /////////////////////////////////////

      if (fabs(pho1->Pt()) < pt1cut)
      {
        continue;
      }

      // another loop to go into the saved cluster
      for (int kCs = 0; kCs < iCs; kCs++)
      {
        if (jCs == kCs)
        {
          continue;
        }

        pho2 = savClusLV[kCs];

        if (fabs(pho2->Pt()) < pt2cut)
        {
          continue;
        }

        alphaCut = fabs((pho1->E() - pho2->E()) / (pho1->E() + pho2->E()));

        if (alphaCut > alphacutval)
        {
          continue;
        }

        TLorentzVector pi0lv;

        if (pho1->DeltaR(*pho2) > deltaRconecut)
        {
          continue;
        }

        pi0lv = *pho1 + *pho2;
        if (fabs(pi0lv.Pt()) > pi0ptcut)
        {
          float pairInvMass = pi0lv.M();

          // fill the tower by tower histograms with invariant mass
          // cemc_hist_eta_phi[_maxTowerEtas[jCs]][_maxTowerPhis[jCs]]->Fill(pairInvMass);
          // not useful in summer 23 data
          eta_hist.at(_maxTowerEtas[jCs])->Fill(pairInvMass);
          pt1_ptpi0_alpha->Fill(pho1->Pt(), pi0lv.Pt(), alphaCut);
          pairInvMassTotal->Fill(pairInvMass);
          mass_eta->Fill(pairInvMass, _clusterEtas[jCs]);
          mass_eta_phi->Fill(pairInvMass, _clusterEtas[jCs], _clusterPhis[jCs]);
        }
      }
    }
  }
  std::cout << "total number of events: " << nEntries << std::endl;
  std::cout << "total number of events discarded: " << discarded_clusters << std::endl;
}

//__________oo00oo__________oo00oo_________________
// This one is for etaslices
void CaloCalibEmc_Pi0::Loop_for_eta_slices(int nevts, const std::string &filename, TTree *intree, const std::string &incorrFile)
{
  //  std::arrays have their indices backward, this is the old float myaggcorr[96][260];
  std::array<std::array<float, 260>, 96> myaggcorr{};
  std::for_each(myaggcorr.begin(), myaggcorr.end(), [](auto &row)
                { row.fill(1.); });

  std::cout << "running w/ corr file? : " << incorrFile << std::endl;

  if (!incorrFile.empty())
  {
    TFile *infileNt = new TFile(incorrFile.c_str());
    std::cout << "loaded incorrFile " << infileNt << std::endl;

    float myieta;
    float myiphi;
    float mycorr;
    float myaggcv;

    TNtuple *innt_corrVals = (TNtuple *) infileNt->Get("nt_corrVals");

    innt_corrVals->SetBranchAddress("tower_eta", &myieta);
    innt_corrVals->SetBranchAddress("tower_phi", &myiphi);
    innt_corrVals->SetBranchAddress("corr_val", &mycorr);
    innt_corrVals->SetBranchAddress("agg_cv", &myaggcv);

    int ntCorrs = innt_corrVals->GetEntries();

    for (int ij = 0; ij < ntCorrs; ij++)
    {
      innt_corrVals->GetEntry(ij);
      int ci = (int) myieta;
      int cj = (int) myiphi;
      myaggcorr.at(ci).at(cj) = myaggcv;
      if (ij > ntCorrs - 2)
      {
        std::cout << "loaded corrs eta,phi,aggcv " << myieta
                  << " " << myiphi << " " << myaggcv << std::endl;
      }
    }

    infileNt->Close();
    delete infileNt;
  }

  std::cout << "in loop" << std::endl;

  TTree *t1 = intree;
  if (!intree)
  {
    TFile *f = new TFile(filename.c_str());
    t1 = (TTree *) f->Get("_eventTree");
  }

  // Set Branches
  //  t1->SetBranchAddress("_eventNumber", &_eventNumber);
  t1->SetBranchAddress("_nClusters", &_nClusters);
  //  t1->SetBranchAddress("_clusterIDs", _clusterIDs);
  t1->SetBranchAddress("_clusterEnergies", _clusterEnergies);
  t1->SetBranchAddress("_clusterPts", _clusterPts);
  t1->SetBranchAddress("_clusterEtas", _clusterEtas);
  t1->SetBranchAddress("_clusterPhis", _clusterPhis);
  t1->SetBranchAddress("_maxTowerEtas", _maxTowerEtas);
  t1->SetBranchAddress("_maxTowerPhis", _maxTowerPhis);

  // pre-loop to save all the clusters LorentzVector

  std::array<TLorentzVector *, 10000> savClusLV{};

  //  int nEntries = (int) t1->GetEntriesFast();
  int nEntries = (int) t1->GetEntries();
  int nevts2 = nevts;

  if (nevts < 0 || nEntries < nevts)
  {
    nevts2 = nEntries;
  }

  for (int i = 0; i < nevts2; i++)
  {
    // load the ith instance of the TTree
    t1->GetEntry(i);

    if ((i % 10 == 0 && i < 200) || (i % 100 == 0 && i < 1000) || (i % 1000 == 0 && i < 37003) || i % 10000 == 0)
    {
      std::cout << "evt no " << i << std::endl;
    }
    // calibration correction will be applied here

    int nClusters = _nClusters;

    if (nClusters > 60)
    {
      continue;
    }

    for (int j = 0; j < nClusters; j++)
    {
      // float px, py, pz;
      float pt, eta, phi, E, aggcv;
      pt = _clusterPts[j];
      eta = _clusterEtas[j];
      phi = _clusterPhis[j];
      E = _clusterEnergies[j];
      aggcv = myaggcorr.at(_maxTowerEtas[j]).at(_maxTowerPhis[j]);

      pt *= aggcv;
      E *= aggcv;

      savClusLV.at(j) = new TLorentzVector();
      savClusLV.at(j)->SetPtEtaPhiE(pt, eta, phi, E);
    }

    TLorentzVector *pho1, *pho2;
    int iCs = nClusters;
    for (int jCs = 0; jCs < iCs; jCs++)
    {
      pho1 = savClusLV.at(jCs);

      if (fabs(pho1->Pt()) < 1.0)
      {
        continue;
      }

      // another loop to go into the saved cluster
      for (int kCs = 0; kCs < iCs; kCs++)
      {
        if (jCs == kCs)
        {
          continue;
        }

        pho2 = savClusLV.at(kCs);

        if (fabs(pho2->Pt()) < 0.6)
        {
          continue;
        }

        TLorentzVector pi0lv;
        if (pho1->DeltaR(*pho2) > 0.45)
        {
          continue;
        }
        pi0lv = *pho1 + *pho2;
        float pairInvMass = pi0lv.M();
        if (pi0lv.Pt() < 1.0)
        {
          continue;
        }

        /*
        if (_maxTowerEtas[jCs]==50 && ((pi0lv.M()>=0.015 && pi0lv.M()<=0.095) || (pi0lv.M()>=0.175 && pi0lv.M()<=0.255)))
        {
                e1_hist_wo_alpha->Fill(pho1->E());
                e2_hist_wo_alpha->Fill(pho2->E());
        }
        */

        alphaCut = fabs((pho1->E() - pho2->E()) / (pho1->E() + pho2->E()));
        if (alphaCut > 0.50)
        {
          continue;  // 0.50 to begin with
        }

        /*
        if  (_maxTowerEtas[jCs]==50 && ((pi0lv.M()>=0.015 && pi0lv.M()<=0.095) || (pi0lv.M()>=0.175 && pi0lv.M()<=0.255)))
        {
                e1_hist_w_alpha->Fill(pho1->E());
                e2_hist_w_alpha->Fill(pho2->E());
        }
        */

        // fill the tower by tower histograms with invariant mass
        // we don't need to fill tower-by-tower level when we do for eta slices
        // although filling here just so we don't have to change codes in other places
        cemc_hist_eta_phi.at(_maxTowerEtas[jCs]).at(_maxTowerPhis[jCs])->Fill(pairInvMass);
        eta_hist.at(_maxTowerEtas[jCs])->Fill(pairInvMass);
        // pt1_ptpi0_alpha->Fill(pho1->Pt(), pi0lv.Pt(), alphaCut);
      }
    }
  }
}

// _______________________________________________________________..
void CaloCalibEmc_Pi0::Fit_Histos(const std::string &incorrFile)
{
  std::cout << " Inside Fit_Histos_Eta_Phi." << std::endl;

  //  std::arrays have their indices backward, this is the old float myaggcorr[96][256];
  std::array<std::array<float, 256>, 96> myaggcorr{};
  std::for_each(myaggcorr.begin(), myaggcorr.end(), [](auto &row)
                { row.fill(1.); });

  if (!incorrFile.empty())
  {
    TFile *infileNt = new TFile(incorrFile.c_str());

    float myieta;
    float myiphi;
    float mycorr;
    float myaggcv;

    TNtuple *innt_corrVals = (TNtuple *) infileNt->Get("nt_corrVals");

    innt_corrVals->SetBranchAddress("tower_eta", &myieta);
    innt_corrVals->SetBranchAddress("tower_phi", &myiphi);
    innt_corrVals->SetBranchAddress("corr_val", &mycorr);
    innt_corrVals->SetBranchAddress("agg_cv", &myaggcv);

    int ntCorrs = innt_corrVals->GetEntries();

    for (int ij = 0; ij < ntCorrs; ij++)
    {
      innt_corrVals->GetEntry(ij);
      int ci = (int) myieta;
      int cj = (int) myiphi;
      myaggcorr.at(ci).at(cj) = myaggcv;
      if (ij > ntCorrs - 2)
      {
        std::cout << "loaded corrs eta,phi,aggcv " << myieta
                  << " " << myiphi << " " << myaggcv << std::endl;
      }
    }

    infileNt->Close();
    delete infileNt;
  }

  cal_output->cd();

  TF1 *f1[25000];
  TF1 *f2[25000];
  TF1 *total[25000];
  TF1 *fit_fn[25000];
  int kj = 0;

  // arrays to hold the fit results (cemc)
  fitp1_eta_phi2d = new TH2F("fitp1_eta_phi2d", "fit p1 eta phi", 96, 0, 96, 256, 0, 256);

  double cemc_par1_values[96][256] = {{0.0}};  // mising braces Werror w/o double braces
  // double cemc_par0_values[96][256] = {0.0};
  // double cemc_par0_errors[96][256] = {0.0};
  double cemc_par1_errors[96][256] = {{0.0}};
  // double cemc_par2_values[96][256] = {0.0};
  // double cemc_par2_errors[96][256] = {0.0};

  // create Ntuple object of the fit result from the data
  TNtuple *nt_corrVals = new TNtuple("nt_corrVals", "Ntuple of the corrections", "tower_eta:tower_phi:corr_val:agg_cv");

  for (int ieta = 0; ieta < 96; ieta++)  // eta loop
  {
    for (int iphi = 0; iphi < 256; iphi++)
    {
      // for(int ithirds=0; ithirds<3; ithirds++)
      //{
      //	for (int ieta=0+ithirds*32; ieta<(ithirds*32+16); ieta++)
      //	{
      //		for (int iphi=0; iphi<16; iphi++)
      //		{
      //		}
      //	}
      // }
      //  find max bin around peak
      float pkloc = 0.0;
      float bsavloc = 0.0;
      for (int kfi = 1; kfi < 20; kfi++)  // old kfi<25
      {
        float locbv = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinContent(kfi);
        if (locbv > bsavloc)
        {
          pkloc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinCenter(kfi);
          bsavloc = locbv;
        }
      }

      f1[kj] = new TF1("f1", "gaus", 0.06, 0.20);  //"gaus",pkloc-0.03,pkloc+0.03
      f2[kj] = new TF1("f2", "pol2", 0.01, 0.4);

      cemc_hist_eta_phi.at(ieta).at(iphi)->Fit(f1[kj], "", "", pkloc - 0.04, pkloc + 0.04);
      float fpkloc2 = f1[kj]->GetParameter(1);

      TGraphErrors *grtemp = new TGraphErrors();
      std::string bkgNm = std::string("grBkgEta_phi_") + std::to_string(ieta) + std::string("_") + std::to_string(iphi);

      std::cout << " getting " << bkgNm << " mean was " << fpkloc2
                << " " << pkloc << std::endl;

      grtemp->SetName(bkgNm.c_str());
      int ingr = 0;
      for (int gj = 1; gj < cemc_hist_eta_phi.at(ieta).at(iphi)->GetNbinsX() + 1; gj++)
      {
        float binc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinCenter(gj);
        float cntc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinContent(gj);
        if ((binc > 0.06 * fpkloc2 / 0.145 && binc < 0.09 * fpkloc2 / 0.145) || (binc > 0.22 * fpkloc2 / 0.145 && binc < 0.35 * fpkloc2 / 0.145))
        {
          grtemp->SetPoint(ingr, binc, cntc);
          grtemp->SetPointError(ingr++, 0.001, sqrt(cntc));
        }
      }

      grtemp->Fit(f2[kj]);

      total[kj] = new TF1("total", "gaus(0)+pol2(3)", 0.06, 0.25);  // 0.3*fpkloc2/0.145

      double par[6];

      f1[kj]->GetParameters(&par[0]);
      f2[kj]->GetParameters(&par[3]);

      total[kj]->SetParameters(par);
      total[kj]->SetParLimits(2, 0.01, 0.027);

      cemc_hist_eta_phi.at(ieta).at(iphi)->Fit(total[kj], "R");
      fit_fn[kj] = cemc_hist_eta_phi.at(ieta).at(iphi)->GetFunction("total");

      grtemp->Write();

      if (fit_fn[kj])
      {
        // cemc_hist_eta_phi.at(ieta).at(iphi) = i;
        cemc_par1_values[ieta][iphi] = fit_fn[kj]->GetParameter(1);

        //	if (!(cemc_par1_values[ieta][iphi]>0.0))
        //	{
        //		cemc_par1_values[ieta][iphi] = 0.5;
        //	}
        // cemc_par0_values[ieta][iphi]	= cemc_eta_phi_result->GetParameter(0);
        cemc_par1_errors[ieta][iphi] = fit_fn[kj]->GetParError(1);
        // cemc_par2_values[ieta][iphi]  = cemc_eta_phi_result->GetParameter(2);
        // cemc_par2_errors[ieta][iphi] = cemc_eta_phi_result->GetParError(2);
      }
      else
      {
        std::cout << "Warning::Fit Failed for eta bin : " << ieta << iphi << std::endl;
        cemc_par1_values[ieta][iphi] = -999.0;
        cemc_par1_errors[ieta][iphi] = -999.0;
      }

      nt_corrVals->Fill(ieta, iphi, 0.135 / cemc_par1_values[ieta][iphi], 0.135 / cemc_par1_values[ieta][iphi] * myaggcorr.at(ieta).at(iphi));
      //}

      // cemc_p1_eta_phi->Fill(cemc_par1_values[ieta][iphi],ieta,iphi);

      // fitp0_eta_phi2d->SetBinContent(ieta+1,iphi+1,cemc_par0_values[ieta][iphi]);
      fitp1_eta_phi2d->SetBinContent(ieta + 1, iphi + 1, cemc_par1_values[ieta][iphi]);
      fitp1_eta_phi2d->SetBinError(ieta + 1, iphi + 1, cemc_par1_errors[ieta][iphi]);

      kj++;
    }
  }

  // nt_corrVals->Fill(ieta,259,0.135/cemc_par1_values[ieta][iphi],0.135/cemc_par1_values[ieta][iphi]*myaggcorr[ieta][259]);

  /*
  TGraphErrors g1(96, eta_value, eta_par1_value, 0, eta_par1_error);
  g1.SetTitle("pi0 mean eta; eta; p1");
  g1.SetMarkerStyle(20);
  g1.SetMarkerColor(2);
  g1.Draw("P");
  g1.SetName("eta_p1");
  g1.Write();

  TGraphErrors g2(96, eta_value, eta_par2_value, 0, eta_par2_error);
  g2.SetTitle("pi0 sigma eta; eta; p2");
  g2.Draw("AP");
  g2.SetName("eta_p2");
  g2.Write();

  */

  fitp1_eta_phi2d->Write();
  std::cout << " finished fit_histos_eta_phi. " << std::endl;
}

// _______________________________________________________________..
void CaloCalibEmc_Pi0::Fit_Histos_Eta_Phi_Add96(const std::string &incorrFile)
{
  std::cout << " Inside Fit_Histos_Eta_Phi." << std::endl;
  //  std::arrays have their indices backward, this is the old float myaggcorr[96][256];
  std::array<std::array<float, 256>, 96> myaggcorr{};
  std::for_each(myaggcorr.begin(), myaggcorr.end(), [](auto &row)
                { row.fill(1.); });

  if (!incorrFile.empty())
  {
    TFile *infileNt = new TFile(incorrFile.c_str());

    float myieta;
    float myiphi;
    float mycorr;
    float myaggcv;

    TNtuple *innt_corrVals = (TNtuple *) infileNt->Get("nt_corrVals");

    innt_corrVals->SetBranchAddress("tower_eta", &myieta);
    innt_corrVals->SetBranchAddress("tower_phi", &myiphi);
    innt_corrVals->SetBranchAddress("corr_val", &mycorr);
    innt_corrVals->SetBranchAddress("agg_cv", &myaggcv);

    int ntCorrs = innt_corrVals->GetEntries();

    for (int ij = 0; ij < ntCorrs; ij++)
    {
      innt_corrVals->GetEntry(ij);
      int ci = (int) myieta;
      int cj = (int) myiphi;
      myaggcorr.at(ci).at(cj) = myaggcv;
      if (ij > ntCorrs - 2)
      {
        std::cout << "loaded corrs eta,phi,aggcv " << myieta
                  << " " << myiphi << " " << myaggcv << std::endl;
      }
    }

    infileNt->Close();
    delete infileNt;
  }

  cal_output->cd();

  TF1 *f1[25000];
  TF1 *f2[25000];
  TF1 *total[25000];
  TF1 *fit_fn[25000];
  int kj = 0;

  // arrays to hold the fit results (cemc)
  fitp1_eta_phi2d = new TH2F("fitp1_eta_phi2d", "fit p1 eta phi", 96, 0, 96, 256, 0, 256);

  double cemc_par1_values[96][256] = {{0.0}};  // mising braces Werror w/o double braces
  // double cemc_par0_values[96][256] = {0.0};
  // double cemc_par0_errors[96][256] = {0.0};
  double cemc_par1_errors[96][256] = {{0.0}};
  // double cemc_par2_values[96][256] = {0.0};
  // double cemc_par2_errors[96][256] = {0.0};

  // create Ntuple object of the fit result from the data
  TNtuple *nt_corrVals = new TNtuple("nt_corrVals", "Ntuple of the corrections", "tower_eta:tower_phi:corr_val:agg_cv");

  for (int ieta = 0; ieta < 96; ieta++)  // eta loop
  {
    for (int iphi = 0; iphi < 256; iphi++)
    {
      if (ieta > 15 || iphi > 15)
      {
        continue;
      }

      // for(int ithirds=0; ithirds<3; ithirds++)
      //{
      //	for (int ieta=0+ithirds*32; ieta<(ithirds*32+16); ieta++)
      //	{
      //		for (int iphi=0; iphi<16; iphi++)
      //		{
      //		}
      //	}
      // }
      //  find max bin around peak
      float pkloc = 0.0;
      float bsavloc = 0.0;
      for (int kfi = 1; kfi < 20; kfi++)  // old kfi<25
      {
        float locbv = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinContent(kfi);
        if (locbv > bsavloc)
        {
          pkloc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinCenter(kfi);
          bsavloc = locbv;
        }
      }

      f1[kj] = new TF1("f1", "gaus", 0.06, 0.20);  //"gaus",pkloc-0.03,pkloc+0.03
      f2[kj] = new TF1("f2", "pol2", 0.01, 0.4);

      cemc_hist_eta_phi.at(ieta).at(iphi)->Fit(f1[kj], "", "", pkloc - 0.04, pkloc + 0.04);
      float fpkloc2 = f1[kj]->GetParameter(1);

      TGraphErrors *grtemp = new TGraphErrors();
      std::string bkgNm = std::string("grBkgEta_phi_") + std::to_string(ieta) + std::string("_") + std::to_string(iphi);

      std::cout << " getting " << bkgNm << " mean was " << fpkloc2
                << " " << pkloc << std::endl;

      grtemp->SetName(bkgNm.c_str());
      int ingr = 0;
      for (int gj = 1; gj < cemc_hist_eta_phi.at(ieta).at(iphi)->GetNbinsX() + 1; gj++)
      {
        float binc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinCenter(gj);
        float cntc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinContent(gj);
        if ((binc > 0.06 * fpkloc2 / 0.145 && binc < 0.09 * fpkloc2 / 0.145) || (binc > 0.22 * fpkloc2 / 0.145 && binc < 0.35 * fpkloc2 / 0.145))
        {
          grtemp->SetPoint(ingr, binc, cntc);
          grtemp->SetPointError(ingr++, 0.001, sqrt(cntc));
        }
      }

      grtemp->Fit(f2[kj]);

      total[kj] = new TF1("total", "gaus(0)+pol2(3)", 0.06, 0.25);  // 0.3*fpkloc2/0.145

      double par[6];

      f1[kj]->GetParameters(&par[0]);
      f2[kj]->GetParameters(&par[3]);

      total[kj]->SetParameters(par);
      total[kj]->SetParLimits(2, 0.01, 0.027);

      cemc_hist_eta_phi.at(ieta).at(iphi)->Fit(total[kj], "R");
      fit_fn[kj] = cemc_hist_eta_phi.at(ieta).at(iphi)->GetFunction("total");

      grtemp->Write();

      if (fit_fn[kj])
      {
        // cemc_hist_eta_phi.at(ieta).at(iphi) = i;
        cemc_par1_values[ieta][iphi] = fit_fn[kj]->GetParameter(1);

        //	if (!(cemc_par1_values[ieta][iphi]>0.0))
        //	{
        //		cemc_par1_values[ieta][iphi] = 0.5;
        //	}
        // cemc_par0_values[ieta][iphi]	= cemc_eta_phi_result->GetParameter(0);
        cemc_par1_errors[ieta][iphi] = fit_fn[kj]->GetParError(1);
        // cemc_par2_values[ieta][iphi]  = cemc_eta_phi_result->GetParameter(2);
        // cemc_par2_errors[ieta][iphi] = cemc_eta_phi_result->GetParError(2);
      }
      else
      {
        std::cout << "Warning::Fit Failed for eta bin : " << ieta << iphi << std::endl;
      }

      for (int ipatt_eta = 0; ipatt_eta < 6; ipatt_eta++)
      {
        for (int ipatt_phi = 0; ipatt_phi < 16; ipatt_phi++)
        {
          // if ((ipatt_eta>0) || (ipatt_phi>0))
          //{
          nt_corrVals->Fill(ieta + ipatt_eta * 16, iphi + ipatt_phi * 16, 0.135 / cemc_par1_values[ieta][iphi], 0.135 / cemc_par1_values[ieta][iphi] * myaggcorr.at(ieta).at(iphi));
          //}
        }
      }

      // nt_corrVals->Fill(ieta,259,0.135/cemc_par1_values[ieta][iphi],0.135/cemc_par1_values[ieta][iphi]*myaggcorr[ieta][259]);

      // cemc_p1_eta_phi->Fill(cemc_par1_values[ieta][iphi],ieta,iphi);

      // fitp0_eta_phi2d->SetBinContent(ieta+1,iphi+1,cemc_par0_values[ieta][iphi]);
      fitp1_eta_phi2d->SetBinContent(ieta + 1, iphi + 1, cemc_par1_values[ieta][iphi]);
      fitp1_eta_phi2d->SetBinError(ieta + 1, iphi + 1, cemc_par1_errors[ieta][iphi]);

      kj++;
    }
  }

  /*
  TGraphErrors g1(96, eta_value, eta_par1_value, 0, eta_par1_error);
  g1.SetTitle("pi0 mean eta; eta; p1");
  g1.SetMarkerStyle(20);
  g1.SetMarkerColor(2);
  g1.Draw("P");
  g1.SetName("eta_p1");
  g1.Write();

  TGraphErrors g2(96, eta_value, eta_par2_value, 0, eta_par2_error);
  g2.SetTitle("pi0 sigma eta; eta; p2");
  g2.Draw("AP");
  g2.SetName("eta_p2");
  g2.Write();

  */

  fitp1_eta_phi2d->Write();
  std::cout << " finished fit_histos_eta_phi. " << std::endl;
}

// _______________________________________________________________..
void CaloCalibEmc_Pi0::Fit_Histos_Eta_Phi_Add32(const std::string &incorrFile)
{
  std::cout << " Inside Fit_Histos_Eta_Phi." << std::endl;

  //  std::arrays have their indices backward, this is the old float myaggcorr[96][256];
  std::array<std::array<float, 256>, 96> myaggcorr{};
  std::for_each(myaggcorr.begin(), myaggcorr.end(), [](auto &row)
                { row.fill(1.); });

  if (!incorrFile.empty())
  {
    TFile *infileNt = new TFile(incorrFile.c_str());

    float myieta;
    float myiphi;
    float mycorr;
    float myaggcv;

    TNtuple *innt_corrVals = (TNtuple *) infileNt->Get("nt_corrVals");

    innt_corrVals->SetBranchAddress("tower_eta", &myieta);
    innt_corrVals->SetBranchAddress("tower_phi", &myiphi);
    innt_corrVals->SetBranchAddress("corr_val", &mycorr);
    innt_corrVals->SetBranchAddress("agg_cv", &myaggcv);

    int ntCorrs = innt_corrVals->GetEntries();

    for (int ij = 0; ij < ntCorrs; ij++)
    {
      innt_corrVals->GetEntry(ij);
      int ci = (int) myieta;
      int cj = (int) myiphi;
      myaggcorr.at(ci).at(cj) = myaggcv;
      if (ij > ntCorrs - 2)
      {
        std::cout << "loaded corrs eta,phi,aggcv " << myieta
                  << " " << myiphi << " " << myaggcv << std::endl;
      }
    }

    infileNt->Close();
    delete infileNt;
  }

  cal_output->cd();

  TF1 *f1[25000];
  TF1 *f2[25000];
  TF1 *total[25000];
  int kj = 0;

  // arrays to hold the fit results (cemc)
  TF1 *cemc_eta_phi_result = nullptr;
  fitp1_eta_phi2d = new TH2F("fitp1_eta_phi2d", "fit p1 eta phi", 96, 0, 96, 256, 0, 256);
  double cemc_par1_values[96][256] = {{0.0}};
  // double cemc_par0_values[96][256] = {0.0};
  // double cemc_par0_errors[96][256] = {0.0};
  double cemc_par1_errors[96][256] = {{0.0}};
  // double cemc_par2_values[96][256] = {0.0};
  // double cemc_par2_errors[96][256] = {0.0};

  // create Ntuple object of the fit result from the data
  TNtuple *nt_corrVals = new TNtuple("nt_corrVals", "Ntuple of the corrections", "tower_eta:tower_phi:corr_val:agg_cv");

  for (int ithirds = 0; ithirds < 3; ithirds++)
  {
    for (int ieta = 0 + ithirds * 32; ieta < (ithirds * 32 + 16); ieta++)
    {
      for (int iphi = 0; iphi < 16; iphi++)
      {
        float pkloc = 0.0;
        float bsavloc = 0.0;
        for (int kfi = 1; kfi < 25; kfi++)  // assuming pi0 peak within 25 bins and no other peak within
        {
          float locbv = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinContent(kfi);
          if (locbv > bsavloc)
          {
            pkloc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinCenter(kfi);
            bsavloc = locbv;
          }
        }

        f1[kj] = new TF1("f1", "gaus", pkloc - 0.03, pkloc + 0.03);
        f2[kj] = new TF1("f2", "pol2", 0.01, 0.4);

        cemc_hist_eta_phi.at(ieta).at(iphi)->Fit(f1[kj], "", "", pkloc - 0.04, pkloc + 0.04);
        float fpkloc2 = f1[kj]->GetParameter(1);

        TGraphErrors *grtemp = new TGraphErrors();
        std::string bkgNm = std::string("grBkgEta_phi_") + std::to_string(ieta) + std::string("_") + std::to_string(iphi);

        std::cout << " getting " << bkgNm << " mean was " << fpkloc2
                  << " " << pkloc << std::endl;

        grtemp->SetName(bkgNm.c_str());
        int ingr = 0;
        for (int gj = 1; gj < cemc_hist_eta_phi.at(ieta).at(iphi)->GetNbinsX() + 1; gj++)
        {
          float binc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinCenter(gj);
          float cntc = cemc_hist_eta_phi.at(ieta).at(iphi)->GetBinContent(gj);
          if ((binc > 0.06 * fpkloc2 / 0.145 && binc < 0.09 * fpkloc2 / 0.145) || (binc > 0.22 * fpkloc2 / 0.145 && binc < 0.35 * fpkloc2 / 0.145))
          {
            grtemp->SetPoint(ingr, binc, cntc);
            grtemp->SetPointError(ingr++, 0.001, sqrt(cntc));
          }
        }

        grtemp->Fit(f2[kj]);

        total[kj] = new TF1("total", "gaus(0)+pol2(3)", 0.06, 0.3 * fpkloc2 / 0.145);

        double par[6];

        f1[kj]->GetParameters(&par[0]);
        f2[kj]->GetParameters(&par[3]);

        total[kj]->SetParameters(par);
        total[kj]->SetParLimits(2, 0.01, 0.027);
        cemc_hist_eta_phi.at(ieta).at(iphi)->Fit(total[kj], "R");
        cemc_eta_phi_result = cemc_hist_eta_phi.at(ieta).at(iphi)->GetFunction("total");
        kj++;

        grtemp->Write();

        if (cemc_eta_phi_result)
        {
          // cemc_hist_eta_phi.at(ieta).at(iphi) = i;
          cemc_par1_values[ieta][iphi] = cemc_eta_phi_result->GetParameter(1);

          //	if (!(cemc_par1_values[ieta][iphi]>0.0))
          //	{
          //		cemc_par1_values[ieta][iphi] = 0.5;
          //	}
          // cemc_par0_values[ieta][iphi]	= cemc_eta_phi_result->GetParameter(0);
          cemc_par1_errors[ieta][iphi] = cemc_eta_phi_result->GetParError(1);
          // cemc_par2_values[ieta][iphi]  = cemc_eta_phi_result->GetParameter(2);
          // cemc_par2_errors[ieta][iphi] = cemc_eta_phi_result->GetParError(2);
        }
        else
        {
          std::cout << "Warning::Fit Failed for eta bin : " << ieta << iphi << std::endl;
        }

        for (int ipatt_eta = 0; ipatt_eta < 2; ipatt_eta++)
        {
          for (int ipatt_phi = 0; ipatt_phi < 16; ipatt_phi++)
          {
            // if ((ipatt_eta>0) || (ipatt_phi>0))
            //{
            nt_corrVals->Fill(ieta + ipatt_eta * 16, iphi + ipatt_phi * 16, 0.135 / cemc_par1_values[ieta][iphi], 0.135 / cemc_par1_values[ieta][iphi] * myaggcorr.at(ieta).at(iphi));
            //}
          }
        }

        // nt_corrVals->Fill(ieta,259,0.135/cemc_par1_values[ieta][iphi],0.135/cemc_par1_values[ieta][iphi]*myaggcorr[ieta][259]);

        // cemc_p1_eta_phi->Fill(cemc_par1_values[ieta][iphi],ieta,iphi);

        // fitp0_eta_phi2d->SetBinContent(ieta+1,iphi+1,cemc_par0_values[ieta][iphi]);
        fitp1_eta_phi2d->SetBinContent(ieta + 1, iphi + 1, cemc_par1_values[ieta][iphi]);
        fitp1_eta_phi2d->SetBinError(ieta + 1, iphi + 1, cemc_par1_errors[ieta][iphi]);
      }
    }
  }
  std::cout << " finished fit_histos_eta_phi. " << std::endl;
}

// _______________________________________________________________..
void CaloCalibEmc_Pi0::Fit_Histos_Etas96(const std::string &incorrFile)
{
  //  std::arrays have their indices backward, this is the old float myaggcorr[96][260];
  std::array<std::array<float, 260>, 96> myaggcorr{};
  std::for_each(myaggcorr.begin(), myaggcorr.end(), [](auto &row)
                { row.fill(1.); });

  if (!incorrFile.empty())
  {
    TFile *infileNt = new TFile(incorrFile.c_str());

    float myieta;
    float myiphi;
    float mycorr;
    float myaggcv;

    TNtuple *innt_corrVals = (TNtuple *) infileNt->Get("nt_corrVals");

    innt_corrVals->SetBranchAddress("tower_eta", &myieta);
    innt_corrVals->SetBranchAddress("tower_phi", &myiphi);
    innt_corrVals->SetBranchAddress("corr_val", &mycorr);
    innt_corrVals->SetBranchAddress("agg_cv", &myaggcv);

    int ntCorrs = innt_corrVals->GetEntries();

    for (int ij = 0; ij < ntCorrs; ij++)
    {
      innt_corrVals->GetEntry(ij);
      int ci = (int) myieta;
      int cj = (int) myiphi;
      myaggcorr.at(ci).at(cj) = myaggcv;
      if (ij > ntCorrs - 2)
      {
        std::cout << "loaded corrs eta,phi,aggcv " << myieta
                  << " " << myiphi << " " << myaggcv << std::endl;
      }
    }

    infileNt->Close();
    delete infileNt;
  }

  cal_output->cd();

  TF1 *f1[96];
  TF1 *f2[96];
  TF1 *total[96];
  TF1 *fit_fn[96];
  int kj = 0;

  // arrays to hold parameters and its error
  double eta_value[96];
  double eta_par1_value[96] = {0.0};
  double eta_par1_error[96] = {0.0};
  double eta_par2_value[96] = {0.0};
  double eta_par2_error[96] = {0.0};
  //  double eta_[96] = {0.0};

  // create Ntuple object from your calculated data
  TNtuple *nt_corrVals = new TNtuple("nt_corrVals", "Ntuple of the corrections", "tower_eta:tower_phi:corr_val:agg_cv");

  float avgmean = 0.170;

  int okHist[96] = {0};

  for (int i = 0; i < 96; i++)
  {
    if (eta_hist.at(i) == nullptr)
    {
      continue;
    }
    else if (eta_hist.at(i)->GetEntries() == 0)
    {
      continue;
    }

    okHist[i] = 1;

    // find max bin around peak
    float pkloc = 0.0;
    float bsavloc = 0.0;

    for (int kfi = 1; kfi < 30; kfi++)
    {
      float locbv = eta_hist.at(i)->GetBinContent(kfi);
      if (locbv > bsavloc)
      {
        pkloc = eta_hist.at(i)->GetBinCenter(kfi);
        bsavloc = locbv;
      }
    }

    f1[kj] = new TF1("f1", "gaus", 0.06, 0.20);  // pkloc-0.03, pkloc+0.03);
    f2[kj] = new TF1("f2", "pol2", 0.005, 0.5);

    if (pkloc < 0.110 || pkloc > 0.2)
    {
      pkloc = 0.170;
    }

    eta_hist.at(i)->Fit(f1[kj], "", "", pkloc - 0.04, pkloc + 0.04);
    float fpkloc2 = f1[kj]->GetParameter(1);

    if (fpkloc2 < 0.110 || fpkloc2 > 0.2)
    {
      fpkloc2 = 0.170;
    }

    TGraphErrors *grtemp = new TGraphErrors();
    std::string bkgNm = std::string("grBkgEta_") + std::to_string(i);

    std::cout << " getting " << bkgNm << "  mean was "
              << fpkloc2 << " " << pkloc << std::endl;
    grtemp->SetName(bkgNm.c_str());
    int ingr = 0;
    int firstnonzbin = -1;
    for (int gj = 1; gj < eta_hist.at(i)->GetNbinsX() + 1; gj++)
    {
      float binc = eta_hist.at(i)->GetBinCenter(gj);
      float cntc = eta_hist.at(i)->GetBinContent(gj);
      if (
          ((binc > 0.03 * fpkloc2 / 0.145 && binc < 0.09 * fpkloc2 / 0.145) || (binc > 0.265 * fpkloc2 / 0.145 && binc < 0.31 * fpkloc2 / 0.145)) && (cntc > 0 || firstnonzbin > 0))
      {
        grtemp->SetPoint(ingr, binc, cntc);
        grtemp->SetPointError(ingr++, 0.001, sqrt(cntc));
      }
      if (cntc && firstnonzbin < 0)
      {
        firstnonzbin = gj;
      }
    }

    grtemp->Fit(f2[kj]);

    total[kj] = new TF1("total", "gaus(0)+pol2(3)", 0.03 * fpkloc2 / 0.145, 0.31 * fpkloc2 / 0.145);  // 0.3*fpkloc2/0.145);

    double par[6];

    f1[kj]->GetParameters(&par[0]);
    f2[kj]->GetParameters(&par[3]);

    total[kj]->SetParameters(par);
    total[kj]->SetParLimits(2, 0.015, 0.037);
    //    total[kj]->SetParLimits(0,0.01*eta_hist.at(i)->Integral(0,40),eta_hist.at(i)->Integral(0,40));
    total[kj]->SetParLimits(0, 5.0, eta_hist.at(i)->Integral(0, 40));
    total[kj]->SetParLimits(1, 0.08, 0.22);
    //    total[kj]->SetParLimits(3,0.01,0.027);
    //    total[kj]->SetParLimits(3,0.01,0.027);

    eta_hist.at(i)->Fit(total[kj], "R");

    fit_fn[kj] = eta_hist.at(i)->GetFunction("total");

    grtemp->Write();

    if (fit_fn[kj])
    {
      eta_value[i] = i;
      eta_par1_value[i] = fit_fn[kj]->GetParameter(1);
      if (!(eta_par1_value[i] > 0))
      {
        eta_par1_value[i] = 0.5;
      }
      eta_par1_error[i] = fit_fn[kj]->GetParError(1);
      eta_par2_value[i] = fit_fn[kj]->GetParameter(2);
      eta_par2_error[i] = fit_fn[kj]->GetParError(2);
    }
    else
    {
      std::cout << "WarninG :: Fit Failed for eta bin : " << i << std::endl;
    }

    avgmean += eta_par1_value[i];

    kj++;
  }

  avgmean /= kj + 1;
  int kj2 = 0;

  for (int iijk = 0; iijk < 96; iijk++)
  {
    if (okHist[iijk] == 0)
    {
      continue;
    }

    float locavg = 0.01;
    if (kj2 < 6 && iijk < 88)
    {
      locavg = (eta_par1_value[iijk + 3] + eta_par1_value[iijk + 4] + eta_par1_value[iijk + 5] + eta_par1_value[iijk + 6] + eta_par1_value[iijk + 7]) / 5.0;
      if (fabs(locavg - avgmean) > 0.2 * avgmean)
      {
        locavg = avgmean;
      }
    }
    else if (iijk > 90)
    {
      locavg = (eta_par1_value[iijk - 3] + eta_par1_value[iijk - 4] + eta_par1_value[iijk - 5] + eta_par1_value[iijk - 6] + eta_par1_value[iijk - 7]) / 5.0;
      if (fabs(locavg - avgmean) > 0.2 * avgmean)
      {
        locavg = avgmean;
      }
    }
    else
    {
      locavg = (eta_par1_value[iijk - 3] + eta_par1_value[iijk - 2] + eta_par1_value[iijk - 1] + eta_par1_value[iijk + 1] + eta_par1_value[iijk + 2] + eta_par1_value[iijk + 3]) / 6.0;
      if (fabs(locavg - avgmean) > 0.2 * avgmean)
      {
        locavg = avgmean;
      }
    }

    float thisfitp1 = 0;
    if (fit_fn[kj2])
    {
      thisfitp1 = fit_fn[kj2]->GetParameter(1);
      if (fabs(thisfitp1 - locavg) > 0.12 * locavg)
      {
        std::cout << "redoing  fit for eta" << iijk << " " << fit_fn[kj2] << std::endl;
        fit_fn[kj2]->SetParameter(1, locavg);
        fit_fn[kj2]->FixParameter(1, locavg);
        eta_hist.at(iijk)->Fit(fit_fn[kj2]);
        fit_fn[kj2]->SetParLimits(1, locavg - 0.011, locavg + 0.011);
        eta_hist.at(iijk)->Fit(fit_fn[kj2]);
      }
      // if (iijk > 63 && iijk < 72)
      //   {
      //     std::cout << "redoing  fit for eta" <<  iijk << " " <<  fit_fn[kj2] << std::endl;
      //     fit_fn[kj2]->SetParameter(1,0.26);
      //     fit_fn[kj2]->FixParameter(1,0.26);
      //     eta_hist.at(iijk)->Fit(fit_fn[kj2]);
      //     //	      fit_fn[kj2]->SetParLimits(1,locavg - 0.011, locavg + 0.011);
      //     //eta_hist.at(iijk)->Fit(fit_fn[kj2]);

      //   }

      eta_value[iijk] = iijk;
      eta_par1_value[iijk] = fit_fn[kj2]->GetParameter(1);
      if (!(eta_par1_value[iijk] > 0))
      {
        eta_par1_value[iijk] = 0.5;
      }
      eta_par1_error[iijk] = fit_fn[kj2]->GetParError(1);
      eta_par2_value[iijk] = fit_fn[kj2]->GetParameter(2);
      eta_par2_error[iijk] = fit_fn[kj2]->GetParError(2);
    }
    else
    {
      std::cout << "WarninG :: Fit Failed for eta bin : " << iijk << std::endl;
    }

    for (int jj = 0; jj < 256; jj++)
    {
      nt_corrVals->Fill(iijk, jj, _setMassVal / eta_par1_value[iijk], _setMassVal / eta_par1_value[iijk] * myaggcorr.at(iijk).at(jj));
    }

    // nt_corrVals->Fill(i,259,0.135/eta_par1_value[i],0.135/eta_par1_value[i]*myaggcorr[i][259]);

    kj2++;

  }  // iijk

  TGraphErrors *g1 = new TGraphErrors(96, eta_value, eta_par1_value, nullptr, eta_par1_error);
  g1->SetTitle("pi0 mean eta; eta; p1");
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(2);
  g1->SetName("eta_p1");

  TGraphErrors *g2 = new TGraphErrors(96, eta_value, eta_par2_value, nullptr, eta_par2_error);
  g2->SetTitle("pi0 sigma eta; eta; p2");
  g2->SetName("eta_p2");

  cal_output->WriteTObject(g1);
  cal_output->WriteTObject(g2);

  std::cout << " finished fit_histos_eta_slice" << std::endl;
}

//_________________________________________________________________________..
void CaloCalibEmc_Pi0::Get_Histos(const std::string &infile, const std::string &outfile)
{
  std::string ts = "cp -rp " + infile + std::string(" ") + outfile;
  gSystem->Exec(ts.c_str());

  cal_output = new TFile(outfile.c_str(), "UPDATE");
  // load the file from the fun4all 1st run

  for (int i = 0; i < 96; i++)
  {
    // getting eta towers
    std::string b = std::string("eta_") + std::to_string(i);
    TH1F *heta_temp = (TH1F *) cal_output->Get(b.c_str());
    eta_hist.at(i) = heta_temp;

    std::cout << "got " << b << std::endl;

    // getting eta_phi towers
    for (int j = 0; j < 256; j++)
    {
      std::string hist_name = std::string("emc_ieta") + std::to_string(i) + std::string("_phi") + std::to_string(j);
      TH1F *h_eta_phi_temp = (TH1F *) cal_output->Get(hist_name.c_str());
      cemc_hist_eta_phi.at(i).at(j) = h_eta_phi_temp;
    }
  }
  std::cout << "finished Loading histos" << std::endl;
}

// _______________________________________________________________________________..
void CaloCalibEmc_Pi0::Add_32()
{
  for (int ithirds = 0; ithirds < 3; ithirds++)
  {
    for (int ieta = 0 + ithirds * 32; ieta < (ithirds * 32 + 16); ieta++)
    {
      for (int iphi = 0; iphi < 16; iphi++)
      {
        for (int ipatt_eta = 0; ipatt_eta < 2; ipatt_eta++)
        {
          for (int ipatt_phi = 0; ipatt_phi < 16; ipatt_phi++)
          {
            if ((ipatt_eta > 0) || (ipatt_phi > 0))
            {
              cemc_hist_eta_phi.at(ieta).at(iphi)->Add(cemc_hist_eta_phi.at(ieta + ipatt_eta * 16).at(iphi + ipatt_phi * 16));
            }
          }
        }
      }
    }
  }
}

// ________________________________________________________________________________..
void CaloCalibEmc_Pi0::Add_96()
{
  std::cout << " Inside Add_96()." << std::endl;
  for (int ieta = 0; ieta < 16; ieta++)
  {
    for (int iphi = 0; iphi < 16; iphi++)
    {
      for (int ipatt_eta = 0; ipatt_eta < 6; ipatt_eta++)
      {
        for (int ipatt_phi = 0; ipatt_phi < 16; ipatt_phi++)
        {
          if ((ipatt_eta > 0) || (ipatt_phi > 0))
          {
            cemc_hist_eta_phi.at(ieta).at(iphi)->Add(cemc_hist_eta_phi.at(ieta + ipatt_eta * 16).at(iphi + ipatt_phi * 16));
          }
        }
      }
    }
  }
  std::cout << " Finished Add_96(). " << std::endl;
}
