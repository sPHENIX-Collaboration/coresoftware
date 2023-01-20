#include "LiteCaloEval.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

class PHCompositeNode;

/// This function is used for the histo fitting process. x is a 1d array that holds xaxis values
/// par is an array of 1d array of parameters we set our fit function to. So p[0] = p[1] = 1 unless 
/// otherwise determined 

// cppcheck suggests to use const which is correct but root
// barfs when this is used with TF1
// cppcheck-suppress constParameter
double LCE_fitf(Double_t *x, Double_t *par)
{
  return par[0] * LCE_grff->Eval(x[0] * par[1], 0, "S");
}

//____________________________________________________________________________..
LiteCaloEval::LiteCaloEval(const std::string &name, const std::string &caloname, const std::string &filename)
  : SubsysReco(name)
  , f_temp(0)
  , fitmin(0.)
  , fitmax(0.)
  , mode(0)
  , _caloname(caloname)
  , _filename(filename)
{
}

///setters
void LiteCaloEval::setFitMax(float fitMax)
{
  fitmax = fitMax;
}

void LiteCaloEval::setFitMin(float fitMin)
{
  fitmin = fitMin;
}


///getters
float LiteCaloEval::getFitMax()
{
  return fitmax;
}

float LiteCaloEval::getFitMin()
{
  return fitmin;
}


//____________________________________________________________________________..
int LiteCaloEval::InitRun(PHCompositeNode * /*topNode*/)
{
  std::cout << "in init run " << std::endl;
  // just quit if we forgot to set the calorimeter type
  if (calotype == LiteCaloEval::NONE)
  {
    std::cout << Name() << ": No Calo Type set" << std::endl;
    gSystem->Exit(1);
  }

  /*
  mygaus = new TF1("mygaus","gaus",-4,4);
  mygaus->SetParameters(1e7,0.,1);
  */

  _ievent = 0;

  cal_output = new TFile(_filename.c_str(), "RECREATE");

  if (calotype == LiteCaloEval::HCALIN)
  {
    hcalin_energy_eta = new TH2F("hcalin_energy_eta", "hcalin energy eta", 1000, 0, 100, 240, -1.1, 1.1);
    hcalin_e_eta_phi = new TH3F("hcalin_e_eta_phi", "hcalin e eta phi", 50, 0, 14, 24, -1.1, 1.1, 64, -M_PI, M_PI);
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_in_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_in_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_in_energy", 40000, 0., 4.);
      }
    }

    for (int i = 0; i < 24; i++)
    {

      

      std::string hist_name = "hcalin_eta_" + std::to_string(i);
      if (i < 24)
	{
	  
	  hcalin_eta[i] = new TH1F(hist_name.c_str(), "hcalin eta's", 40000, 0, 4.);
	}
      else
	{
	  hcalin_eta[i] = new TH1F(hist_name.c_str(), "hcalin eta's", 2000000, 0, 4.);
	}
    }

  }
  else if (calotype == LiteCaloEval::HCALOUT)
  {
    hcalout_energy_eta = new TH2F("hcalout_energy_eta", "hcalout energy eta", 100, 0, 10, 240, -1.1, 1.1);
    hcalout_e_eta_phi = new TH3F("hcalout_e_eta_phi", "hcalout e eta phi", 48, 0, 10, 24, -1.1, 1.1, 64, -M_PI, M_PI);
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_out_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_out_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_out energy", 50000, 0, 10.0);
      }
    }

    for (int i = 0; i < 25; i++)
    {
      std::string hist_name = "hcalout_eta_" + std::to_string(i);
      if (i < 24)
      {
        hcalout_eta[i] = new TH1F(hist_name.c_str(), "hcalout eta's", 50000, 0, 10.);
      }
      else
      {
        hcalout_eta[i] = new TH1F(hist_name.c_str(), "hcalout eta's", 1000000, 0, 10);
      }
    }
  }
  else if (calotype == LiteCaloEval::CEMC)
  {
    for (int i = 0; i < 96; i++)
    {
      for (int j = 0; j < 258; j++)
      {
        std::string hist_name = "emc_ieta" + std::to_string(i) + "_phi" + std::to_string(j);

        cemc_hist_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hist_ieta_phi_leaf(e)", 4000, 0, 8.0);
      }
    }

    //create the eta 1d histos
    for (int i = 0; i < 97; i++)
    {
      gStyle->SetOptFit(1);
      std::string b = "eta_" + std::to_string(i);

      if (i < 96)
      {
        eta_hist[i] = new TH1F(b.c_str(), "eta and all phi's", 4000, 0, 8.0);
      }
      else
      {
        eta_hist[i] = new TH1F(b.c_str(), "eta and all phi's", 1000000, 0, 10);
      }
    }

    //make 2d histo
    energy_eta_hist = new TH2F("energy_eta_hist", "energy eta and all phi", 512, 0, 10, 960, -1.15, 1.15);

    //make 3d histo
    e_eta_phi = new TH3F("e_eta_phi", "e v eta v phi", 50, 0, 10, 192, -1.1335, 1.13350, 256, -M_PI, M_PI);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::process_event(PHCompositeNode *topNode)
{
  if (_ievent % 100 == 0) std::cout << "LiteCaloEval::process_event(PHCompositeNode *topNode) Processing Event " << _ievent << std::endl;

  std::string towernode = "TOWER_CALIB_" + _caloname;
  RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode, towernode);
  if (!towers)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
    exit(-1);
  }

  std::string towergeomnode = "TOWERGEOM_" + _caloname;
  RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode);
  if (!towergeom)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << towergeomnode << std::endl;
    exit(-1);
  }

  if (mode && _ievent < 4)
  {
    std::cout << "mode is set " << std::endl;
    //   e*= 1.15;
  }

  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;

    if (tower->get_energy() < 0.0005) continue;

    RawTowerGeom *tower_geom = towergeom->get_tower_geometry(tower->get_id());
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
    float e = tower->get_energy();

    if (calotype == LiteCaloEval::CEMC)
    {
      if (towerid < 0)
      {
        std::cout << "a towerid was less than 0 " << std::endl;
        break;
      }

      //      e = e + 0.3*mygaus->GetRandom();
      // if (mode)
      // 	{
      // 	  e*= 1.15;
      // 	}

      //fill the hist with energy data
      cemc_hist_eta_phi[ieta][iphi]->Fill(e);

      //fill and fit the 1d hist eta and all phi
      eta_hist[96]->Fill(e);
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

      // if (mode)
      // 	{
      // 	  //e = e + 0.27*mygaus->GetRandom()*fabs(1-0.12/e);
      // 	  int ket = ieta/2;
      // 	  //     int ll = (6-abs(6-k));
      // 	  int llet = ket%6;
      // 	  e *= 0.945+llet*0.02;
      // 	  int pket = iphi/4;
      // 	  if (pket%2==0)
      // 	    e*= 1.03;
      // 	}

      hcal_out_eta_phi[ieta][iphi]->Fill(e);

      hcalout_eta[ieta]->Fill(e);
      hcalout_eta[24]->Fill(e);

      hcalout_energy_eta->Fill(e, eta);

      hcalout_e_eta_phi->Fill(e, eta, phi);
    }
    else if (calotype == LiteCaloEval::HCALIN)
    {
      //fill the hist with energy data

      hcal_in_eta_phi[ieta][iphi]->Fill(e);

      hcalin_eta[ieta]->Fill(e);
      hcalin_eta[24]->Fill(e);

      hcalin_energy_eta->Fill(e, eta);

      hcalin_e_eta_phi->Fill(e, eta, phi);
    }
  }

  _ievent++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::End(PHCompositeNode * /*topNode*/)
{
  cal_output->cd();
  std::cout << " writing lite calo file" << std::endl;
  cal_output->Write();
  //cout <<" wrote lite calo file" << endl;

  return Fun4AllReturnCodes::EVENT_OK;

  /* 
  // turn off fitting at end because
     // it always needs merged first
  
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
  else if (calotype == LiteCaloEval::CEMC && 0)
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

  std::cout <<" writing lite calo file" << std::endl;
  cal_output->Write();
  //cout <<" wrote lite calo file" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
  */
}

void LiteCaloEval::Get_Histos(const char *infile, const char *outfile)
{
  std::string outF = outfile;
  std::string inF = infile;
  if ((inF == ""))
  {
    std::cout << "need infile to run LiteCaloEval::Get_Histos" << std::endl;
    exit(0);
  }

  if (!(outF == ""))
  {
    TString ts = "cp -rp ";
    ts += infile;
    ts += " ";
    ts += outfile;
    gSystem->Exec(ts.Data());
    f_temp = new TFile(outfile, "UPDATE");  // load the file from the fun4all 1st run
  }
  else
  {
    f_temp = new TFile(infile);
  }

  int max_ieta = 96;
  if (calotype != LiteCaloEval::CEMC)
  {
    max_ieta = 24;
  }

  int max_iphi = 256;
  if (calotype != LiteCaloEval::CEMC)
  {
    max_iphi = 64;
  }

  for (int i = 0; i < max_ieta; i++)
  {
    TString a;
    a.Form("%d", i);
    TString b = "eta_" + a;
    if (calotype == LiteCaloEval::HCALOUT)
    {
      b = "hcalout_" + b;
    }
    else if(calotype == LiteCaloEval::HCALIN)
      b = "hcalin_" + b;

      /// holds the eta slice of histos
    TH1F *heta_temp = (TH1F *) f_temp->Get(b.Data());
    if (!heta_temp && i == 0) std::cout << " warning, hist " << b.Data() << " not found" << std::endl;

    /// assign the eta slice histo to an array (these arrays are private members in LCE.h)
    eta_hist[i] = heta_temp;
    if (calotype == LiteCaloEval::HCALOUT)
      hcalout_eta[i] = heta_temp;
    else if (calotype == LiteCaloEval::HCALIN)
      hcalin_eta[i] = heta_temp;

    for (int j = 0; j < max_iphi; j++)
    {
      std::string hist_name_p = "emc_ieta" + std::to_string(i) + "_phi" + std::to_string(j);
      if (calotype == LiteCaloEval::HCALOUT)
        hist_name_p = "hcal_out_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);
      else if (calotype == LiteCaloEval::HCALIN)
	hist_name_p = "hcal_in_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

      if (j < 2)
      {
        std::cout << "getting " << hist_name_p << std::endl;
      }
      TH1F *heta_tempp = (TH1F *) f_temp->Get(hist_name_p.c_str());

      if (!heta_tempp && i == 0)
      {
        std::cout << " warning hist " << hist_name_p.c_str() << " not found" << std::endl;
      }

      cemc_hist_eta_phi[i][j] = heta_tempp;
      if (calotype == LiteCaloEval::HCALOUT)
        hcal_out_eta_phi[i][j] = heta_tempp;
      else if (calotype == LiteCaloEval::HCALIN)
	hcal_in_eta_phi[i][j] = heta_tempp;
      
    }
  }
}

void LiteCaloEval::FitRelativeShifts(LiteCaloEval *ref_lce, int modeFitShifts)
{
  bool onlyEta = false;
  if (fitmin < 0.001)
  {
    fitmin = 0.15;
  }
  if (fitmax < 0.001)
  {
    fitmax = 1.3;
  }
  float par_value[96];
  float par_err[96];
  float eta_value[96];
  float eta_err[96];
  //  float rel_err[96];

  if (f_temp) f_temp->cd();

  TF1 *f1 = new TF1("myexpo", LCE_fitf, 0.1, 10, 2);
  //  TF1 * myexpo = new TF1("myexpo","expo",0.1,10);

  f1->SetParameters(1e2, 1.0);

  if (modeFitShifts % 10 == 1)
  {
    onlyEta = true;
  }

  /// nsmooth is an automatic smoothing of histos
  int nsmooth = 1;

/// if want more smoothing, look at mFS. If mFS=10, then addSmooth=1, and then there will be 2 total smoothings 
  int addSmooth = modeFitShifts % 100 / 10;
  if (addSmooth)
  {
    nsmooth += addSmooth;
  }

  /// flag that will run fits at tower level to get shifts, or will run at ring level and fit towers to eta slice histo
  bool flag_fit_rings =false;

  /// e.x. if mFS = 010 , 010%1000 = 10, 10/100 = 0  -> false, so dont run tower to rings
  /// e.x. if mFS = 110, 110%1000 = 110, 110/100 = 1 -> true, so run to fit eta slice histos
  if (modeFitShifts % 1000 / 100 == 1)
    flag_fit_rings = true;

  int max_ieta = 96;
  if (calotype != LiteCaloEval::CEMC)
  {
    max_ieta = 24;
  }

  int max_iphi = 256;
  if (calotype != LiteCaloEval::CEMC)
  {
    max_iphi = 64;
  }

  TH2F *corrPat = new TH2F("corrPat", "", max_ieta, 0, max_ieta, max_iphi, 0, max_iphi);

  int minbin = 0, maxbin = max_ieta;

  //  if (calotype == LiteCaloEval::CEMC)
  //{
  //  minbin = 94;  //maxbin = 49;
  // }

  for (int i = minbin; i < maxbin; i++)
  {
    //a,b hold the p1 value, p1 error, and then repeats for the second fit (c,d) and third fit (e,f)

    double a1f;
    double b1f;

    TString myClnm = "newhc_eta";
    myClnm = myClnm + i;
    TString myClnm2 = "gnewhc_eta";
    myClnm2 = myClnm2 + i;
    std::cout << " making fitting  " << myClnm.Data() << std::endl;

    int iik = i;
    TH1F *hnewf = 0;
    if (calotype == LiteCaloEval::CEMC)
      hnewf = (TH1F *) ref_lce->eta_hist[iik]->Clone(myClnm.Data());
    
    else if (calotype == LiteCaloEval::HCALOUT)
      hnewf = (TH1F *) ref_lce->hcalout_eta[iik]->Clone(myClnm.Data());
    
    else if (calotype == LiteCaloEval::HCALIN)
      hnewf = (TH1F *) ref_lce->hcalin_eta[iik]->Clone(myClnm.Data()); 
    
    std::cout << "got neweff " << std::endl;
    hnewf->Smooth(1);
    if (nsmooth > 1) 
      hnewf->Smooth(nsmooth-1);
    
    LCE_grff = new TGraph(hnewf);

    f1->SetParameters(1e1, 1.0);

    TF1 *f2f = 0;
    if (calotype == LiteCaloEval::CEMC)
    {
      if (nsmooth > 1)
      {
        eta_hist[i]->Smooth(nsmooth);
      }
      eta_hist[i]->Fit("myexpo", "L", "", fitmin, fitmax);
      //hcalout_eta[i]->Fit("f2","R+");
      //hcalout_eta[i]->Fit("f3","R+");

      f2f = (TF1 *) eta_hist[i]->GetFunction("myexpo");
    }

    else if (calotype == LiteCaloEval::HCALOUT)
    {
      if (i == 0) std::cout << "in hcal" << std::endl;
      if (nsmooth > 1)
      {
        hcalout_eta[i]->Smooth(nsmooth);
      }
      hcalout_eta[i]->Fit("myexpo", "L", "", fitmin, fitmax);
      //hcalout_eta[i]->Fit("f2","R+");
      //hcalout_eta[i]->Fit("f3","R+");

      f2f = (TF1 *) hcalout_eta[i]->GetFunction("myexpo");
    }

    else if (calotype == LiteCaloEval::HCALIN)
      {
	//if (i == 0) std::cout << "in hcalin" << std::endl;
	if (nsmooth > 1)
	  hcalin_eta[i]->Smooth(nsmooth);
	
	hcalin_eta[i]->Fit("myexpo","L","",fitmin,fitmax);
	
	f2f = (TF1 *) hcalin_eta[i]->GetFunction("myexpo");
      }
    

    TGraph *grT = new TGraph(1000);
    grT->SetName(myClnm2.Data());
    for (int jjk = 0; jjk < 1000; jjk++)
    {
      float xjj = fitmin + jjk * (fitmax - fitmin) / 1000.0;
      grT->SetPoint(jjk, xjj, f2f->Eval(xjj));
    }

    grT->Write();

    a1f = f2f->GetParameter(1);
    b1f = f2f->GetParError(1);

    //assign retreived parameter values
    par_value[i] = a1f;
    par_err[i] = b1f;
    eta_value[i] = i;
    eta_err[i] = 0.01;
    //       rel_err[i] = par_err[i] / par_value[i];


    if (onlyEta) continue;

    // else also continue with each individual phi tower

    for (int j = 0; j < max_iphi; j++)
    {
      TString myClnmp = "newhc_eta";
      myClnmp = myClnmp + 1000 * (i + 2) + j;
      TString myClnm2p = "gnewhc_eta";
      myClnm2p = myClnm2p + 1000 * (i + 2) + j;
      if (j == 0)
      {
        std::cout << " making fitting  " << myClnmp.Data() << std::endl;
      }
      //       int iikp = i;
      //  int jjkp = i;
 
      /// histo to hold tower from modified root files
      TH1F * hnewfp  = 0;
      /// histo to hold eta slice with current tower removed
      TH1F * rmTower = 0;
      
      if (calotype == LiteCaloEval::CEMC)
      {
        hnewfp = (TH1F *) ref_lce->cemc_hist_eta_phi[i][j]->Clone(myClnmp.Data());
	if(flag_fit_rings == true)
	  {
	    rmTower = (TH1F *)ref_lce->eta_hist[i]->Clone();
	    rmTower->Add((TH1F *)ref_lce->cemc_hist_eta_phi[i][j],-1.0);
	  }
      }
      else if (calotype == LiteCaloEval::HCALOUT)
	{
        hnewfp = (TH1F *) ref_lce->hcal_out_eta_phi[i][j]->Clone(myClnmp.Data());
	if(flag_fit_rings == true)
	  {
	    rmTower = (TH1F *)ref_lce->hcalout_eta[i]->Clone();
	    rmTower->Add((TH1F *)ref_lce->hcal_out_eta_phi[i][j],-1.0);
	  }
	}
      else if (calotype == LiteCaloEval::HCALIN)
	{
	  hnewfp = (TH1F *) ref_lce->hcal_in_eta_phi[i][j]->Clone(myClnmp.Data());
	  
	  if(flag_fit_rings == true)
	    {
	      rmTower = (TH1F *)ref_lce->hcalin_eta[i]->Clone();
	      rmTower->Add((TH1F *)ref_lce->hcal_in_eta_phi[i][j],-1.0);
	    }
	}
      
      if (j < 3)
      {
        std::cout << "got neweff phi eta  ... fitting w/ fit min,max: " << fitmin
                  << " " << fitmax << std::endl;
      }

      hnewfp->Smooth(nsmooth);
      /*
       hnewfp->Rebin(4);
       hnewfp->Smooth(nsmooth);
       */
      //       hnewfp->Rebin(2);

      //       hnewfp->Smooth(1);

      /*
	If false make tgraph out of tower.
	
	If true (when we fit the towers) they are actually fit against the eta ring. This happends bc "myexpo" uses a tgraph 
	to fit, and if we dont make LCE_grff for a tower, myexpo uses the latest version of LCE_grff, which is an eta slice tgraph
      */

      LCE_grff = new TGraph(hnewfp);
	  
      if(flag_fit_rings == true)
	LCE_grff = new TGraph(rmTower);
      

      f1->SetParameters(1.0, 1.0);

      TF1 *f2f2 = 0;
      if (calotype == LiteCaloEval::CEMC)
      {
        //	   if (nsmooth > 1 )
        //  cemc_hist_eta_phi[i][j]->Smooth(nsmooth);
        // cemc_hist_eta_phi[i][j]->Rebin(2);
        //	   cemc_hist_eta_phi[i][j]->Smooth(1);
        if (j < 2)
	  {
	    cemc_hist_eta_phi[i][j]->Fit("myexpo", "L", "", fitmin, fitmax);
	  }
        else
	  {
	    cemc_hist_eta_phi[i][j]->Fit("myexpo", "LQ", "", fitmin, fitmax);
	  }
        //hcalout_eta[i]->Fit("f2","R+");
        //hcalout_eta[i]->Fit("f3","R+");

        f2f2 = (TF1 *) cemc_hist_eta_phi[i][j]->GetFunction("myexpo");
      }
      else if (calotype == LiteCaloEval::HCALOUT)
      {
        if (i == 0) std::cout << "in hcal" << std::endl;
        hcal_out_eta_phi[i][j]->Smooth(nsmooth);
        hcal_out_eta_phi[i][j]->Rebin(4);
        hcal_out_eta_phi[i][j]->Smooth(nsmooth);

        if (j < 2)
        {
          hcal_out_eta_phi[i][j]->Fit("myexpo", "L", "", fitmin, fitmax);
        }
        else
        {
          hcal_out_eta_phi[i][j]->Fit("myexpo", "LQ", "", fitmin, fitmax);
        }
        //hcalout_eta[i]->Fit("f2","R+");
        //hcalout_eta[i]->Fit("f3","R+");

        f2f2 = (TF1 *) hcal_out_eta_phi[i][j]->GetFunction("myexpo");
      }
      else if (calotype == LiteCaloEval::HCALIN)
      {
	     
	if (j < 2)
	  {
	    hcal_in_eta_phi[i][j]->Fit("myexpo","L","",fitmin,fitmax);
	  }
	else
	  {
	    hcal_in_eta_phi[i][j]->Fit("myexpo","LQ","",fitmin,fitmax);
	  }
	
	f2f2 = (TF1 *)hcal_in_eta_phi[i][j]->GetFunction("myexpo");
	
      }
      

      TGraph *grT2 = new TGraph(1000);
      grT2->SetName(myClnm2p.Data());
      for (int jjk = 0; jjk < 1000; jjk++)
      {
        float xjj = fitmin + jjk * (fitmax - fitmin) / 1000.0;
        grT2->SetPoint(jjk, xjj, f2f->Eval(xjj));
      }

      grT2->Write();

      float a1fp = f2f2->GetParameter(1);
      float b1fp = f2f->GetParError(1);
      corrPat->SetBinContent(i + 1, j + 1, a1fp);
      corrPat->SetBinError(i + 1, j + 1, b1fp);
    } // phi loop
  }//i eta loop

  //create graphs
  TGraphErrors g1(max_ieta, eta_value, par_value, eta_err, par_err);
  g1.SetTitle("fitted shifts; eta; p1");
  g1.SetMarkerStyle(20);
  g1.Draw("ap");
  g1.SetName("Fit1_etaout");
  g1.Write();

  corrPat->Write();

  if (calotype == LiteCaloEval::CEMC)
  {
    std::cout << "TowerSlope module:  writing emcal correction tree into output file"
              << std::endl;
    TTree *t1 = new TTree("emc_corr_tree", "a tree of simple emcal calib corrections");
    int towid;
    float corr;
    t1->Branch("corr", &corr, "corr/F");
    t1->Branch("towid", &towid, "towid/I");

    for (int mjl = 0; mjl < max_ieta; mjl++)
    {
      for (int mjk = 0; mjk < max_iphi; mjk++)
      {
        towid = mjl * 1000 + mjk;
        corr = corrPat->GetBinContent(mjl + 1, mjk + 1);
        if (!(corr > 0.0))
        {
          corr = 1.0;
        }
        else
        {
          corr = 1.0 / corr;
        }
        t1->Fill();
      }
    }
    t1->Write();
  }

  if (calotype == LiteCaloEval::HCALOUT ||
      calotype == LiteCaloEval::HCALIN)
  {
    std::string hcal_corr_file_name = "HCAL_CORR_TXTFILE";

    if (f_temp)
    {
      hcal_corr_file_name += f_temp->GetName();
      if (calotype == LiteCaloEval::HCALOUT)
	hcal_corr_file_name += "OUT" ;
      else
	hcal_corr_file_name += "IN";
      
      hcal_corr_file_name += ".txt";
    }

    std::cout << "TowerSlope module:  writing hcal corrections into output file "
              << hcal_corr_file_name
              << std::endl;
    
    std::ofstream out_hcal_corrF(hcal_corr_file_name.c_str());

    for (int mjl = 0; mjl < max_ieta; mjl++)
      {
	for (int mjk = 0; mjk < max_iphi; mjk++)
	  {
        float corr = corrPat->GetBinContent(mjl + 1, mjk + 1);
        if (!(corr > 0.))
	  {
          corr = 1.0;
        }
        else
        {
          corr = 1.0 / corr;
        }
        out_hcal_corrF << mjl << " "
                       << mjk << " "
                       << corr << std::endl;
      }
    }

    out_hcal_corrF.close();
  }

  if (f_temp)
  {
    f_temp->Write();
    f_temp->Close();
  }
}
