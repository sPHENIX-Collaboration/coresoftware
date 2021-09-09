#include "LiteCaloEval.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <memory>

//____________________________________________________________________________..
LiteCaloEval::LiteCaloEval(const std::string& name, const std::string& caloname, const std::string& filename)
  : SubsysReco(name)
  , _caloname(caloname)
  , _filename(filename)
{
}

//____________________________________________________________________________..
int LiteCaloEval::InitRun(PHCompositeNode* /*topNode*/)
{
  // just quit if we forgot to set the calorimeter type
  if (calotype == LiteCaloEval::NONE)
  {
    std::cout << Name() << ": No Calo Type set" << std::endl;
    gSystem->Exit(1);
  }

  _ievent = 0;

  cal_output = new TFile(_filename.c_str(), "RECREATE");

  if (calotype == LiteCaloEval::HCALIN)
  {
    hcalin_energy_eta = new TH2F("hcalin_energy_eta", "hcalin energy eta", 1000, 0, 10, 240, -1.1, 1.1);
    hcalin_e_eta_phi = new TH3F("hcalin_e_eta_phi", "hcalin e eta phi", 50, 0, 10, 24, -1.1, 1.1, 64, -3.14159, 3.14159);
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_in_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_in_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_in_energy", 1000, 0, 10);
      }
    }

    for (int i = 0; i < 24; i++)
    {
      std::string hist_name = "hcalin_eta_" + std::to_string(i);

      hcalin_eta[i] = new TH1F(hist_name.c_str(), "hcalin eta's", 1000, 0, 10);
    }
  }
  else if (calotype == LiteCaloEval::HCALOUT)
  {
    hcalout_energy_eta = new TH2F("hcalout_energy_eta", "hcalout energy eta", 10, 0, 10, 24000, -1.1, 1.1);
    hcalout_e_eta_phi = new TH3F("hcalout_e_eta_phi", "hcalout e eta phi", 50, 0, 10, 24, -1.1, 1.1, 64, -3.14159, 3.14159);
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_out_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_out_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_out energy", 1000, 0, 10);
      }
    }

    for (int i = 0; i < 24; i++)
    {
      std::string hist_name = "hcalout_eta_" + std::to_string(i);

      hcalout_eta[i] = new TH1F(hist_name.c_str(), "hcalout eta's", 1000, 0, 10);
    }
  }
  else if (calotype == LiteCaloEval::CEMC)
  {
    for (int i = 0; i < 96; i++)
    {
      for (int j = 0; j < 258; j++)
      {
        std::string hist_name = "emc_ieta" + std::to_string(i) + "_phi" + std::to_string(j);

        cemc_hist_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hist_ieta_phi_leaf(e)", 1000, 0, 10);
      }
    }

    //create the eta 1d histos
    for (int i = 0; i < 96; i++)
    {
      gStyle->SetOptFit(1);
      std::string b = "eta_" + std::to_string(i);

      eta_hist[i] = new TH1F(b.c_str(), "eta and all phi's", 1000, 0, 10);
    }

    //make 2d histo
    energy_eta_hist = new TH2F("energy_eta_hist", "energy eta and all phi", 10, 0, 10, 9600, -1, 1);

    //make 3d histo
    e_eta_phi = new TH3F("e_eta_phi", "e v eta v phi", 50, 0, 10, 100, -1, 1, 256, -3.14159, 3.14159);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::process_event(PHCompositeNode* topNode)
{
  if (_ievent % 100 == 0) std::cout << "LiteCaloEval::process_event(PHCompositeNode *topNode) Processing Event " << _ievent << std::endl;

  std::string towernode = "TOWER_CALIB_" + _caloname;
  RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, towernode);
  if (!towers)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
    exit(-1);
  }

  std::string towergeomnode = "TOWERGEOM_" + _caloname;
  RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode);
  if (!towergeom)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << towergeomnode << std::endl;
    exit(-1);
  }

  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    RawTower* tower = rtiter->second;

    if (tower->get_energy() < 0.0005) continue;

    RawTowerGeom* tower_geom = towergeom->get_tower_geometry(tower->get_id());
    if (!tower_geom)
    {
      std::cout << PHWHERE << " ERROR: Can't find tower geometry for this tower hit: ";
      tower->identify();
      exit(-1);
    }

    const int towerid = tower->get_id();
    const int ieta = tower->get_bineta();
    const int iphi = tower->get_binphi();
    const float eta = tower_geom->get_eta();
    const float phi = tower_geom->get_phi();
    const float e = tower->get_energy();
    if (calotype == LiteCaloEval::CEMC)
    {
      if (towerid < 0)
      {
        std::cout << "a towerid was less than 0 " << std::endl;
        break;
      }

      //fill the hist with energy data
      cemc_hist_eta_phi[ieta][iphi]->Fill(e);

      //fill and fit the 1d hist eta and all phi
      eta_hist[ieta]->Fill(e);

      //fill the 2d histo eta, energy and all phi
      energy_eta_hist->Fill(e, eta);

      //fill 3d histo e_eta_phid
      e_eta_phi->Fill(e, eta, phi);
    }
    else if (calotype == LiteCaloEval::HCALOUT)
    {
      //fill the hist with energy data
      //std::cout << ieta << " " <<  iphi  << std::endl;

      hcal_out_eta_phi[ieta][iphi]->Fill(e);

      hcalout_eta[ieta]->Fill(e);

      hcalout_energy_eta->Fill(e, eta);

      hcalout_e_eta_phi->Fill(e, eta, phi);
    }
    else if (calotype == LiteCaloEval::HCALIN)
    {
      //fill the hist with energy data

      hcal_in_eta_phi[ieta][iphi]->Fill(e);

      hcalin_eta[ieta]->Fill(e);

      hcalin_energy_eta->Fill(e, eta);

      hcalin_e_eta_phi->Fill(e, eta, phi);
    }
  }

  _ievent++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::End(PHCompositeNode* /*topNode*/)
{
  cal_output->cd();
  if (calotype == LiteCaloEval::HCALIN)
  {
    double par_value[24];
    double eta_value[24];
    double par_err[24];
    double rel_err[24];

    //same as above but for the second fit performed below in code
    double par_value2[24];
    double par_err2[24];
    double rel_err2[24];

    for (int i = 0; i < 24; i++)
    {
      //a,b,c,d hold the par value, error, par2 value, par2 error respectivley. The 2 refers to the second fit below
      double a;
      double b;
      double c;
      double d;

      //make functions to fit the eta histos
      std::unique_ptr<TF1> f1(new TF1("f1", "expo", 0.02, 0.1));
      std::unique_ptr<TF1> f2(new TF1("f2", "expo", 0.1, 1));
      f2->SetLineColor(7);

      hcalin_eta[i]->Fit("f1", "R");
      hcalin_eta[i]->Fit("f2", "R+");

      a = f1->GetParameter(1);
      b = f1->GetParError(1);

      c = f2->GetParameter(1);
      d = f2->GetParError(1);

      //assign retreived parameter values
      par_value[i] = abs(a);
      par_err[i] = b;
      eta_value[i] = i;
      rel_err[i] = par_err[i] / par_value[i];

      par_value2[i] = abs(c);
      par_err2[i] = d;
      rel_err2[i] = par_err2[i] / par_value2[i];
    }

    //create graphs
    TGraph g1(24, eta_value, par_value);
    g1.SetTitle("HCal In (0.02-0.1 GeV); eta; p1");
    g1.SetMarkerStyle(20);
    g1.Draw("ap");
    g1.SetName("Fit1_hcalin");
    g1.Write();

    TGraph g2(24, eta_value, rel_err);
    g2.SetTitle("HCal In Error (0.02-0.1 GeV); eta; p1 rel error");
    g2.SetMarkerStyle(20);
    g2.Draw("ap");
    g2.SetName("Fit1_err_hcalin");
    g2.Write();

    TGraph g3(24, eta_value, par_value2);
    g3.SetTitle("HCal In (0.1-1 GeV); eta; p1");
    g3.SetMarkerStyle(20);
    g3.Draw("ap");
    g3.SetName("Fit2_hcalin");
    g3.Write();

    TGraph g4(24, eta_value, rel_err2);
    g4.SetTitle("HCal In Error (0.1-1 GeV); eta; p1 rel error");
    g4.SetMarkerStyle(20);
    g4.Draw("ap");
    g4.SetName("Fit2_err_hcalin");
    g4.Write();
  }
  else if (calotype == LiteCaloEval::HCALOUT)
  {
    double par_value[24];
    double eta_value[24];
    double par_err[24];
    double rel_err[24];

    //same as above but for the second/third fit performed below in code
    double par_value2[24];
    double par_err2[24];
    double rel_err2[24];

    double par_value3[24];
    double par_err3[24];
    double rel_err3[24];

    for (int i = 0; i < 24; i++)
    {
      //a,b hold the p1 value, p1 error, and then repeats for the second fit (c,d) and third fit (e,f)
      double a;
      double b;
      double c;
      double d;
      double e;
      double f;

      //make functions to fit the eta histos
      std::unique_ptr<TF1> f1(new TF1("f1", "expo", 0.05, 0.2));
      std::unique_ptr<TF1> f2(new TF1("f2", "expo", 0.2, 1));
      std::unique_ptr<TF1> f3(new TF1("f3", "expo", 1, 2));
      f2->SetLineColor(7);
      f3->SetLineColor(1);

      hcalout_eta[i]->Fit("f1", "R");
      hcalout_eta[i]->Fit("f2", "R+");
      hcalout_eta[i]->Fit("f3", "R+");

      a = f1->GetParameter(1);
      b = f1->GetParError(1);

      c = f2->GetParameter(1);
      d = f2->GetParError(1);

      e = f3->GetParameter(1);
      f = f3->GetParError(1);

      //assign retreived parameter values
      par_value[i] = abs(a);
      par_err[i] = b;
      eta_value[i] = i;
      rel_err[i] = par_err[i] / par_value[i];

      par_value2[i] = abs(c);
      par_err2[i] = d;
      rel_err2[i] = par_err2[i] / par_value2[i];

      par_value3[i] = abs(e);
      par_err3[i] = f;
      rel_err3[i] = par_err3[i] / par_value3[i];
    }

    //create graphs
    TGraph g1(24, eta_value, par_value);
    g1.SetTitle("HCal Out (0.05-0.2 GeV); eta; p1");
    g1.SetMarkerStyle(20);
    g1.Draw("ap");
    g1.SetName("Fit1_hcalout");
    g1.Write();

    TGraph g2(24, eta_value, rel_err);
    g2.SetTitle("HCal Out Error (0.05-0.2 GeV); eta; p1 rel err");
    g2.SetMarkerStyle(20);
    g2.Draw("ap");
    g2.SetName("Fit1_err_hcalout");
    g2.Write();

    TGraph g3(24, eta_value, par_value2);
    g3.SetTitle("HCal Out (0.2-1 GeV); eta; p1");
    g3.SetMarkerStyle(20);
    g3.Draw("ap");
    g3.SetName("Fit2_hcalout");
    g3.Write();

    TGraph g4(24, eta_value, rel_err2);
    g4.SetTitle("HCal Out Error (0.2-1 GeV); eta; p1 rel err");
    g4.SetMarkerStyle(20);
    g4.Draw("ap");
    g4.SetName("Fit2_err_hcalout");
    g4.Write();

    TGraph g5(24, eta_value, par_value3);
    g5.SetTitle("HCal Out (1-2 GeV); eta; p1");
    g5.SetMarkerStyle(20);
    g5.Draw("ap");
    g5.SetName("Fit3_hcalout");
    g5.Write();

    TGraph g6(24, eta_value, rel_err3);
    g6.SetTitle("HCal Out Error (1-2 GeV); eta; p1 rel err");
    g6.SetMarkerStyle(20);
    g6.Draw("ap");
    g6.SetName("Fit3_err_hcalout");
    g6.Write();
  }
  else if (calotype == LiteCaloEval::CEMC)
  {
    //create arrays that holds parameter (p1) value and error
    double par_value[96];
    double eta_value[96];
    double par_err[96];
    double rel_err[96];

    //same as above but for the second fit performed below in code
    double par_value2[96];
    double par_err2[96];
    double rel_err2[96];

    double par_value3[96];
    double par_err3[96];
    double rel_err3[96];

    //create graphs for parameter (p1) vs eta
    for (int i = 0; i < 96; i++)
    {
      //a,b,c,d hold the par value, error, par2 value, par2 error respectivley. The 2 refers to the second fit below
      double a;
      double b;
      double c;
      double d;
      double e;
      double f;

      //make functions to fit the eta histos
      std::unique_ptr<TF1> f1(new TF1("f1", "expo", 0.04, 0.1));
      std::unique_ptr<TF1> f2(new TF1("f2", "expo", 0.1, 0.4));
      std::unique_ptr<TF1> f3(new TF1("f3", "expo", 0.4, 2));
      f2->SetLineColor(7);
      f3->SetLineColor(1);

      eta_hist[i]->Fit("f1", "R");
      eta_hist[i]->Fit("f2", "R+");
      eta_hist[i]->Fit("f3", "R+");

      a = f1->GetParameter(1);
      b = f1->GetParError(1);

      c = f2->GetParameter(1);
      d = f2->GetParError(1);

      e = f3->GetParameter(1);
      f = f3->GetParError(1);

      //assign retreived parameter values
      par_value[i] = abs(a);
      par_err[i] = b;
      eta_value[i] = i;
      rel_err[i] = par_err[i] / par_value[i];

      par_value2[i] = abs(c);
      par_err2[i] = d;
      rel_err2[i] = par_err2[i] / par_value2[i];

      par_value3[i] = abs(e);
      par_err3[i] = f;
      rel_err3[i] = par_err3[i] / par_value3[i];
    }

    //create graphs
    TGraph g1(96, eta_value, par_value);
    g1.SetTitle("EMCal (0.04-0.1 GeV); eta; p1");
    g1.SetMarkerStyle(20);
    g1.Draw("ap");
    g1.SetName("Fit1_emc");
    g1.Write();

    TGraph g2(96, eta_value, rel_err);
    g2.SetTitle("EMCal Error (0.04-0.1 GeV); eta; p1 rel error");
    g2.SetMarkerStyle(20);
    g2.Draw("ap");
    g2.SetName("Fit1_err_emc");
    g2.Write();

    TGraph g3(96, eta_value, par_value2);
    g3.SetTitle("EMCal (0.1-0.4 GeV); eta; p1");
    g3.SetMarkerStyle(20);
    g3.Draw("ap");
    g3.SetName("Fit2_emc");
    g3.Write();

    TGraph g4(96, eta_value, rel_err2);
    g4.SetTitle("EMCal Error (0.1-0.4 GeV); eta; p1 rel error");
    g4.SetMarkerStyle(20);
    g4.Draw("ap");
    g4.SetName("Fit2_err_emc");
    g4.Write();

    TGraph g5(96, eta_value, par_value3);
    g5.SetTitle("EMCal (0.4-2 GeV); eta; p1");
    g5.SetMarkerStyle(20);
    g5.Draw("ap");
    g5.SetName("Fit3_emc");
    g5.Write();

    TGraph g6(96, eta_value, rel_err3);
    g6.SetTitle("EMCal Error (0.4-2 GeV); eta; p1 rel err");
    g6.SetMarkerStyle(20);
    g6.Draw("ap");
    g6.SetName("Fit3_err_emc");
    g6.Write();
  }

  cal_output->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}
