#include "pi0EtaByEta.h"

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// Tower includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TTree.h>

#include <iostream>
#include <sstream>
#include <string>
#include <utility>

pi0EtaByEta::pi0EtaByEta(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("HCALIN")
  , outfilename(filename)
{
  h_mass_eta_lt.fill(nullptr);
  clusMix = new std::vector<std::vector<std::vector<CLHEP::Hep3Vector>>>();
}

pi0EtaByEta::~pi0EtaByEta()
{
  delete hm;
  delete g4hitntuple;
  delete g4cellntuple;
  delete towerntuple;
  delete clusterntuple;
}

int pi0EtaByEta::Init(PHCompositeNode* /*unused*/)
{
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here

  outfile = new TFile(outfilename.c_str(), "RECREATE");

  // correlation plots
  for (int i = 0; i < 96; i++)
  {
    std::string histoname = "h_mass_eta_lt" + std::to_string(i);
    h_mass_eta_lt[i] = new TH1F(histoname.c_str(), "", 50, 0, 0.5);
  }

  h_cemc_etaphi = new TH2F("h_cemc_etaphi", "", 96, 0, 96, 256, 0, 256);
  h_cemc_etaphi_noCalib = new TH2F("h_cemc_etaphi_noCalib", "", 96, 0, 96, 256, 0, 256);

  // 1D distributions
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 120, 0, 1.2);
  h_InvMassMix = new TH1F("h_InvMassMix", "Invariant Mass", 120, 0, 1.2);

  // cluster QA
  h_etaphi_clus = new TH2F("h_etaphi_clus", "", 140, -1.2, 1.2, 64, -1 * M_PI, M_PI);
  h_clus_pt = new TH1F("h_clus_pt", "", 100, 0, 10);

  h_emcal_e_eta = new TH1F("h_emcal_e_eta", "", 96, 0, 96);

  h_pt1 = new TH1F("h_pt1", "", 100, 0, 5);
  h_pt2 = new TH1F("h_pt2", "", 100, 0, 5);

  h_nclusters = new TH1F("h_nclusters", "", 1000, 0, 1000);

  h_nclus_bin = new TH1F("h_nclus_bin", "", NBinsClus, 0, 300);
  h_vtx_bin = new TH1F("h_vtx_bin", "", NBinsVtx, -30, 30);

  for (int ib = 0; ib < NBinsClus; ib++)
  {
    std::string histoname = "h_InvMass_Nclus" + std::to_string(ib);
    h_InvMass_Nclus[ib] = new TH1F(histoname.c_str(), "", 120, 0, 1.2);
  }

  h_event = new TH1F("h_event", "", 1, 0, 1);

  h_pipT_Nclus_mass = new TH3F("h_pipT_Nclus_mass", "", 20, 0, 10, 20, 0, 400, 25, 0, 0.5);

  std::vector<std::vector<CLHEP::Hep3Vector>> temp2 = std::vector<std::vector<CLHEP::Hep3Vector>>();
  std::vector<CLHEP::Hep3Vector> temp = std::vector<CLHEP::Hep3Vector>();
  for (int i = 0; i < NBinsVtx; i++)
  {
    temp2.push_back(temp);
  }
  for (int i = 0; i < NBinsClus; i++)
  {
    clusMix->push_back(temp2);
  }

  return 0;
}

int pi0EtaByEta::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;

  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int pi0EtaByEta::process_towers(PHCompositeNode* topNode)
{
  if ((_eventcounter % 1000) == 0)
  {
    std::cout << _eventcounter << std::endl;
  }

  // float emcaldownscale = 1000000 / 800;

  float emcal_hit_threshold = 0.4;  // GeV

  // cuts
  float maxDr = 1.1;
  float maxAlpha = 0.6;
  float clus_chisq_cut = 4;
  float nClus_ptCut = 0.5;
  int max_nClusCount = 300;

  //----------------------------------get vertex------------------------------------------------------//
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    // std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing"<< std::endl;
    std::cout << "pi0EtaByEta GlobalVertexMap node is missing" << std::endl;
    // return Fun4AllReturnCodes::ABORTRUN;
  }
  float vtx_z = 0;
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex* vtx = vertexmap->begin()->second;
    if (vtx)
    {
      vtx_z = vtx->get_z();
    }
  }

  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (towers)
  {
    int size = towers->size();
    for (int channel = 0; channel < size; channel++)
    {
      TowerInfo* tower = towers->get_tower_at_channel(channel);
      float offlineenergy = tower->get_energy();
      unsigned int towerkey = towers->encode_key(channel);
      int ieta = towers->getTowerEtaBin(towerkey);
      int iphi = towers->getTowerPhiBin(towerkey);
      bool isGood = tower->get_isGood();
      if (offlineenergy > emcal_hit_threshold && isGood)
      {
        h_cemc_etaphi->Fill(ieta, iphi);
      }
      if (tower->get_isNoCalib())
      {
        h_cemc_etaphi_noCalib->Fill(ieta, iphi);
      }
    }
  }

  std::string cluster_node_name = "CLUSTERINFO_CEMC2";
  if (use_pdc)
  {
    cluster_node_name = "CLUSTERINFO_POS_COR_CEMC";
  }
  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, cluster_node_name);
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
    throw std::runtime_error("failed to find TOWERGEOM node in RawClusterDeadHotMask::CreateNodeTree");
  }

  /////////////////////////////
  // clusters
  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;
  RawClusterContainer::ConstIterator clusterIter2;
  int nClusCount = 0;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

    float clus_pt = E_vec_cluster.perp();
    float clus_chisq = recoCluster->get_chi2();

    if (clus_pt < nClus_ptCut)
    {
      continue;
    }
    if (clus_chisq > clus_chisq_cut)
    {
      continue;
    }
    nClusCount++;
  }

  h_nclusters->Fill(nClusCount);

  int nClusBin = h_nclus_bin->FindBin(nClusCount) - 1;
  if (nClusBin == -1 || nClusBin == h_nclus_bin->GetNbinsX())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  int vtxBin = h_vtx_bin->FindBin(vtx_z) - 1;
  if (vtxBin == -1 || vtxBin == h_vtx_bin->GetNbinsX())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  if (vtx_z == -0)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (nClusCount > max_nClusCount)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  float ptClusMax = 7;
  float pt1ClusCut = pt1BaseClusCut;  // 1.3
  float pt2ClusCut = pt2BaseClusCut;  // 0.7

  if (nClusCount > 30)
  {
    pt1ClusCut += NclusDeptFac * (nClusCount - 29) / 200.0;
    pt2ClusCut += NclusDeptFac * (nClusCount - 29) / 200.0;
  }

  float pi0ptcut = 1.22 * (pt1ClusCut + pt2ClusCut);

  h_event->Fill(0);

  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

    float clusE = E_vec_cluster.mag();
    float clus_eta = E_vec_cluster.pseudoRapidity();
    float clus_phi = E_vec_cluster.phi();
    float clus_pt = E_vec_cluster.perp();
    float clus_chisq = recoCluster->get_chi2();

    if (clus_chisq > clus_chisq_cut)
    {
      continue;
    }
    h_clus_pt->Fill(clus_pt);

    // loop over the towers in the cluster
    RawCluster::TowerConstRange towerCR = recoCluster->get_towers();
    RawCluster::TowerConstIterator toweriter;
    float lt_e = -1000;
    unsigned int lt_eta = -1;
    for (toweriter = towerCR.first; toweriter != towerCR.second; ++toweriter)
    {
      int towereta = m_geometry->get_tower_geometry(toweriter->first)->get_bineta();
      int towerphi = m_geometry->get_tower_geometry(toweriter->first)->get_binphi();
      unsigned int key = TowerInfoDefs::encode_emcal(towereta, towerphi);
      unsigned int channel = towers->decode_key(key);
      float energy = towers->get_tower_at_channel(channel)->get_energy();
      if (energy > lt_e)
      {
        lt_e = energy;
        lt_eta = towereta;
      }
    }

    h_etaphi_clus->Fill(clus_eta, clus_phi);

    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);

    if (clus_pt < pt1ClusCut || clus_pt > ptClusMax)
    {
      continue;
    }

    for (clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; clusterIter2++)
    {
      if (clusterIter2 == clusterIter)
      {
        continue;
      }

      RawCluster* recoCluster2 = clusterIter2->second;

      CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*recoCluster2, vertex);

      float clus2E = E_vec_cluster2.mag();
      float clus2_eta = E_vec_cluster2.pseudoRapidity();
      float clus2_phi = E_vec_cluster2.phi();
      float clus2_pt = E_vec_cluster2.perp();
      float clus2_chisq = recoCluster2->get_chi2();

      if (clus2_pt < pt2ClusCut || clus_pt > ptClusMax)
      {
        continue;
      }
      if (clus2_chisq > clus_chisq_cut)
      {
        continue;
      }

      TLorentzVector photon2;
      photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);

      if (fabs(clusE - clus2E) / (clusE + clus2E) > maxAlpha)
      {
        continue;
      }

      if (photon1.DeltaR(photon2) > maxDr)
      {
        continue;
      }

      TLorentzVector pi0 = photon1 + photon2;
      if (pi0.Pt() < pi0ptcut)
      {
        continue;
      }

      h_pt1->Fill(photon1.Pt());
      h_pt2->Fill(photon2.Pt());

      h_InvMass->Fill(pi0.M());
      if (clus2_pt < pt1ClusCut)
      {
        h_InvMass->Fill(pi0.M());
      }
      h_pipT_Nclus_mass->Fill(pi0.Pt(), nClusCount, pi0.M());
      if (lt_eta > 95)
      {
        continue;
      }
      h_mass_eta_lt[lt_eta]->Fill(pi0.M());
      h_InvMass_Nclus[nClusBin]->Fill(pi0.M());
    }

    if (doMix)
    {
      for (const auto& E_vec_cluster2 : clusMix->at(nClusBin)[vtxBin])
      {
        float clus2E = E_vec_cluster2.mag();
        float clus2_eta = E_vec_cluster2.pseudoRapidity();
        float clus2_phi = E_vec_cluster2.phi();
        float clus2_pt = E_vec_cluster2.perp();

        TLorentzVector photon2;
        photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);

        if (fabs(clusE - clus2E) / (clusE + clus2E) > maxAlpha)
        {
          continue;
        }

        if (photon1.DeltaR(photon2) > maxDr)
        {
          continue;
        }

        TLorentzVector pi0 = photon1 + photon2;
        if (pi0.Pt() < pi0ptcut)
        {
          continue;
        }

        h_InvMassMix->Fill(pi0.M());
      }
    }  // doMix
  }    // clus1 loop

  if (doMix)
  {
    clusMix->at(nClusBin)[vtxBin].clear();

    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
    {
      RawCluster* recoCluster = clusterIter->second;

      CLHEP::Hep3Vector vertex(0, 0, vtx_z);
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

      float clus_pt = E_vec_cluster.perp();

      if (recoCluster->get_chi2() > clus_chisq_cut)
      {
        continue;
      }
      if (clus_pt < pt2ClusCut || clus_pt > ptClusMax)
      {
        continue;
      }

      clusMix->at(nClusBin)[vtxBin].push_back(E_vec_cluster);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int pi0EtaByEta::End(PHCompositeNode* /*topNode*/)
{
  outfile->cd();

  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  return 0;
}

TF1* pi0EtaByEta::fitHistogram(TH1* h)
{
  TF1* f_sig_initial = new TF1("f_sig_initial", "[0]/[2]/2.5*exp(-0.5*((x-[1])/[2])^2)", 0.05, 0.25);

  f_sig_initial->SetParLimits(0, 0, 10);
  f_sig_initial->SetParLimits(1, 0.113, 0.2);
  f_sig_initial->SetParLimits(2, 0.007, 0.04);
  f_sig_initial->SetParameter(0, 0.1);
  f_sig_initial->SetParameter(1, 0.15);
  f_sig_initial->SetParameter(2, 0.015);

  h->Fit(f_sig_initial, "N");

  TF1* fitFunc = new TF1("fitFunc", "[0]/[2]/2.5*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2 + [6]*x^3", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0, f_sig_initial->GetParameter(0));
  fitFunc->SetParameter(1, f_sig_initial->GetParameter(1));
  fitFunc->SetParameter(2, f_sig_initial->GetParameter(2));
  fitFunc->SetParameter(3, 0.0);
  fitFunc->SetParameter(4, 0.0);
  fitFunc->SetParameter(5, 0.0);
  fitFunc->SetParameter(6, 0);

  fitFunc->SetParLimits(0, 0, 10);
  fitFunc->SetParLimits(1, 0.113, 0.25);
  fitFunc->SetParLimits(2, 0.007, 0.06);
  fitFunc->SetParLimits(3, -2, 1);
  fitFunc->SetParLimits(4, -1, 40);
  fitFunc->SetParLimits(5, -150, 50);
  // fitFunc->SetParLimits(6, -1,200 );
  fitFunc->SetParLimits(6, -1, 1);

  fitFunc->SetRange(0.05, 0.7);

  // Perform the fit
  h->Fit("fitFunc", "QN");

  return fitFunc;
}

void pi0EtaByEta::fitEtaSlices(const std::string& infile, const std::string& fitOutFile, const std::string& cdbFile)
{
  TFile* fin = new TFile(infile.c_str());
  std::cout << "getting hists" << std::endl;
  TH1F* h_peak_eta = new TH1F("h_peak_eta", "", 96, 0, 96);
  TH1F* h_sigma_eta = new TH1F("h_sigma_eta", "", 96, 0, 96);
  TH1F* h_p3_eta = new TH1F("h_p3_eta", "", 96, 0, 96);
  TH1F* h_p4_eta = new TH1F("h_p4_eta", "", 96, 0, 96);
  TH1F* h_p5_eta = new TH1F("h_p5_eta", "", 96, 0, 96);
  TH1F* h_p6_eta = new TH1F("h_p6_eta", "", 96, 0, 96);
  TH1F* h_p0_eta = new TH1F("h_p0_eta", "", 96, 0, 96);
  if (!fin)
  {
    std::cout << "pi0EtaByEta::fitEtaSlices null fin" << std::endl;
    exit(1);
  }
  TH1F* h_M_eta[96];
  for (int i = 0; i < 96; i++)
  {
    std::string histoname = "h_mass_eta_lt" + std::to_string(i);
    h_M_eta[i] = (TH1F*) fin->Get(histoname.c_str());
    h_M_eta[i]->Scale(1. / h_M_eta[i]->Integral(), "width");
  }

  TF1* fitFunOut[96];
  for (int i = 0; i < 96; i++)
  {
    if (!h_M_eta[i])
    {
      std::cout << "pi0EtaByEta::fitEtaSlices null hist" << std::endl;
    }

    fitFunOut[i] = fitHistogram(h_M_eta[i]);
    std::string funcname = "f_pi0_eta" + std::to_string(i);
    fitFunOut[i]->SetName(funcname.c_str());
    float mass_val_out = fitFunOut[i]->GetParameter(1);
    float mass_err_out = fitFunOut[i]->GetParError(1);
    h_peak_eta->SetBinContent(i + 1, mass_val_out);
    if (isnan(h_M_eta[i]->GetEntries()))
    {
      h_peak_eta->SetBinError(i + 1, 0);
      continue;
    }
    h_peak_eta->SetBinError(i + 1, mass_err_out);
    h_sigma_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(2));
    h_sigma_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(2));
    h_p3_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(3));
    h_p3_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(3));
    h_p4_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(4));
    h_p4_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(4));
    h_p5_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(5));
    h_p5_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(5));
    h_p6_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(6));
    h_p6_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(6));
    h_p0_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(0));
    h_p0_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(0));
  }

  CDBTTree* cdbttree1 = new CDBTTree(cdbFile.c_str());
  CDBTTree* cdbttree2 = new CDBTTree(cdbFile.c_str());

  std::string m_fieldname = "Femc_datadriven_qm1_correction";

  float final_mass_target = target_pi0_mass;

  for (int i = 0; i < 96; i++)
  {
    for (int j = 0; j < 256; j++)
    {
      if (use_h_target_mass)
      {
        final_mass_target = h_target_mass->GetBinContent(i + 1);
      }
      float correction = final_mass_target / h_peak_eta->GetBinContent(i + 1);
      unsigned int key = TowerInfoDefs::encode_emcal(i, j);
      float val1 = cdbttree1->GetFloatValue(key, m_fieldname);
      cdbttree2->SetFloatValue(key, m_fieldname, val1 * correction);
    }
  }

  cdbttree2->Commit();
  cdbttree2->WriteCDBTTree();
  delete cdbttree2;
  delete cdbttree1;

  TFile* fit_out = new TFile(fitOutFile.c_str(), "recreate");
  fit_out->cd();
  for (auto& i : h_M_eta)
  {
    i->Write();
    delete i;
  }
  for (auto& i : fitFunOut)
  {
    i->Write();
    delete i;
  }

  h_p3_eta->Write();
  h_p4_eta->Write();
  h_p5_eta->Write();
  h_p6_eta->Write();
  h_p0_eta->Write();
  h_sigma_eta->Write();
  h_peak_eta->Write();
  fin->Close();

  std::cout << "finish fitting suc" << std::endl;

  return;
}

void pi0EtaByEta::set_massTargetHistFile(const std::string& file)
{
  TFile* fin = new TFile(file.c_str());
  if (fin)
  {
    h_target_mass = (TH1F*) fin->Get("h_massTargetHist");
  }
  else
  {
    std::cout << "pi0EtaByEta::set_massTargetHistFile  No mass file found" << std::endl;
    return;
  }
  if (h_target_mass)
  {
    use_h_target_mass = true;
    std::cout << "using target mass histogram" << std::endl;
  }
  else
  {
    std::cout << "pi0EtaByEta::set_massTargetHistFile  No mass hist found" << std::endl;
  }
  return;
}
