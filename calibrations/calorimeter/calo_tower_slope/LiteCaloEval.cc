#include "LiteCaloEval.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <RtypesCore.h>  // for Double_t
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

class RawTowerGeom;

/*
#include <centrality/CentralityInfo.h>    // for CentralityInfo, CentralityI...
#include <centrality/CentralityInfov1.h>
*/

/// This function is used for the histo fitting process. x is a 1d array that holds xaxis values.
/// par is an array of 1d array of parameters we set our fit function to. So p[0] = p[1] = 1 unless
/// otherwise determined
double LCE_fitf(Double_t *x, Double_t *par)
{
  return par[0] * LCE_grff->Eval(x[0] * par[1], nullptr, "S");
}

//____________________________________________________________________________..
LiteCaloEval::LiteCaloEval(const std::string &name, const std::string &caloname, const std::string &filename)
  : SubsysReco(name)
  , _caloname(caloname)
  , _filename(filename)
  , _inputnodename("TOWERINFO")
  , m_UseTowerInfo(1)
{
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
    hcalin_e_eta_phi = new TH3F("hcalin_e_eta_phi", "hcalin e eta phi", 50, 0, 14, 24, -1.1, 1.1, 64, -3.14159, 3.14159);

    /// create tower histos
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_in_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_in_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_in_energy", 40000, 0, 4);
      }
    }

    for (int i = 0; i < 25; i++)
    {
      std::string hist_name = "hcalin_eta_" + std::to_string(i);
      //      hcalin_eta[i] = new TH1F(hist_name.c_str(), "hcalin eta's", 40000, 0, 4.);
      // now fixed < 25 it should'v been // cppcheck flags this as always true (which it is in this loop
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
    hcalout_e_eta_phi = new TH3F("hcalout_e_eta_phi", "hcalout e eta phi", 48, 0, 10, 24, -1.1, 1.1, 64, -3.14159, 3.14159);

    /// create tower histos
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_out_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_out_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_out energy", 20000, 0, 10);
      }
    }

    /// create eta slice histos
    for (int i = 0; i < 25; i++)
    {
      std::string hist_name = "hcalout_eta_" + std::to_string(i);
      if (i < 24)
      {
        hcalout_eta[i] = new TH1F(hist_name.c_str(), "hcalout eta's", 20000, 0, 10);
      }
      else
      {
        hcalout_eta[i] = new TH1F(hist_name.c_str(), "hcalout eta's", 1000000, 0, 10);
      }
    }
  }

  else if (calotype == LiteCaloEval::CEMC)
  {
    /// create tower histos
    for (int i = 0; i < 96; i++)
    {
      for (int j = 0; j < 256; j++)
      {
        std::string hist_name = "emc_ieta" + std::to_string(i) + "_phi" + std::to_string(j);

        cemc_hist_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hist_ieta_phi_leaf(e)", 5000, 0, 10);
      }
    }

    // create eta slice histos
    for (int i = 0; i < 97; i++)
    {
      gStyle->SetOptFit(1);
      std::string b = "eta_" + std::to_string(i);

      if (i < 96)
      {
        eta_hist[i] = new TH1F(b.c_str(), "eta and all phi's", 5000, 0, 10);
      }
      else
      {
        eta_hist[i] = new TH1F(b.c_str(), "eta and all phi's", 1000000, 0, 10);
      }
    }

    // make 2d histo
    energy_eta_hist = new TH2F("energy_eta_hist", "energy eta and all phi", 512, 0, 10, 960, -1.15, 1.15);

    // make 3d histo
    e_eta_phi = new TH3F("e_eta_phi", "e v eta v phi", 50, 0, 10, 192, -1.1335, 1.13350, 256, -3.14159, 3.14159);
  }

  /*
  //centrality histo
  if(calotype == LiteCaloEval::CEMC)
    evtCentrality = new TH2F("centVsNtow","",100,0,100,24576,0,24576);
  else
    evtCentrality = new TH2F("centVsNtow","",100,0,100,1536,0,1536);
  */

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::process_event(PHCompositeNode *topNode)
{
  if (_ievent % 100 == 0)
  {
    std::cout << "LiteCaloEval::process_event(PHCompositeNode *topNode) Processing Event " << _ievent << std::endl;
  }

  // used to get centrality of collision
  //  we don't need centrality right now
  //  previous code for this was for tests
  /*
  CentralityInfo *cent = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent)
    {
      std::cout << " ERROR -- can't find CentralityInfo node after it should have been created" << std::endl;
      exit(-1);

    }
  */

  // raw tower container
  std::string towernode = "TOWER_CALIB_" + _caloname;
  RawTowerContainer *towers = nullptr;
  RawTowerGeomContainer *towergeom = nullptr;

  // get tower energy for towerslope method only
  RawTowerContainer::ConstRange begin_end;
  RawTowerContainer::ConstIterator rtiter;

  if (m_UseTowerInfo < 1)
  {
    towers = findNode::getClass<RawTowerContainer>(topNode, towernode);
    if (!towers)
    {
      std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
      exit(-1);
    }

    // raw tower geom container
    std::string towergeomnode = "TOWERGEOM_" + _caloname;
    towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode);
    if (!towergeom)
    {
      std::cout << PHWHERE << " ERROR: Can't find " << towergeomnode << std::endl;
      exit(-1);
    }

    begin_end = towers->getTowers();
    rtiter = begin_end.first;
  }

  if (mode && _ievent < 4)
  {
    std::cout << "mode is set " << std::endl;
    //   e*= 1.15;
  }

  //  int centval1 = cent->get_centile(CentralityInfo::PROP::mbd_NS);
  TowerInfoContainer *towerinfos = nullptr;

  // if using towerinfo create a tower object
  if (m_UseTowerInfo)
  {
    //    towernode = "TOWERINFO_CALIB_" + _caloname;
    //    towernode = "TOWERS_Calib_" + _caloname;
    towernode = _inputnodename;

    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, towernode.c_str());

    if (!towerinfos)
    {
      std::cout << PHWHERE << " ERROR: Can't find " << towernode << std::endl;
      exit(-1);
    }
  }

  // getting size of node
  unsigned int nchannels = -1;

  if (m_UseTowerInfo)
  {
    nchannels = towerinfos->size();
  }
  else
  {
    nchannels = towers->size();
  }

  TowerInfo *tower_info;
  RawTower *tower;
  RawTowerGeom *tower_geom;

  //  int towCnt = 0;

  /*
  if (_ievent < 4)
    {
      std::cout << Name() << "towerinfo " << towerinfos << " towers "
                << towers << " nchannels" << nchannels << std::endl;
    }
  */

  for (unsigned int channel = 0; channel < nchannels; channel++)
  {
    if (!m_UseTowerInfo && rtiter == begin_end.second)
    {
      std::cout << "ERROR: Recheck iteration process with rtiter variable" << std::endl;
    }

    float e = -1.0;  // tower energy
    unsigned int ieta = 99999;
    unsigned int iphi = 99999;
    //      int towerid =-1;

    if (m_UseTowerInfo)
    {
      tower_info = towerinfos->get_tower_at_channel(channel);
      unsigned int towerkey = towerinfos->encode_key(channel);

      e = tower_info->get_energy();
      ieta = towerinfos->getTowerEtaBin(towerkey);
      iphi = towerinfos->getTowerPhiBin(towerkey);
      if (!tower_info->get_isGood())
      {
        continue;
      }
    }

    else
    {
      tower = rtiter->second;

      /*
      if (tower->get_energy() < 0.02)
        {
          continue;
        }
      */
      tower_geom = towergeom->get_tower_geometry(tower->get_id());

      if (!tower_geom)
      {
        std::cout << PHWHERE << " ERROR: Can't find tower geometry for this tower hit: ";
        tower->identify();
        exit(-1);
      }

      e = tower->get_energy();
      ieta = tower->get_bineta();
      iphi = tower->get_binphi();

    }  // end else for rawtower mode

    if (ieta > 95 || iphi > 256)
    {
      // rough check for all calos uing the largest emcal
      std::cout << "ieta or iphi not set correctly/ was less than 0 " << std::endl;
      break;
    }

    if (e <= 0)
    {
      continue;
    }

    if (calotype == LiteCaloEval::CEMC)
    {
      /*
      if (e > 0)
        towCnt++;

      if (_ievent > 2 && _ievent < 7 && (channel < 10 || nchannels-channel < 10 || (e > 0.2 && towCnt < nchannels/4.)))
    {
      std::cout << _ievent << " evtn,iphi,ieta, e"  << iphi << " "
                << ieta << " " << e << "  " << Name() << std::endl;
    }
      */

      if (mode)
      {
        int ket = ieta / 8;
        int llet = ket % 6;
        int pket = iphi / 8;
        int ppkket = pket % 3;

        e *= 0.88 + llet * 0.04 - 0.01 + 0.01 * ppkket;
      }

      // fill the hist with energy data
      cemc_hist_eta_phi[ieta][iphi]->Fill(e);

      // fill the 1d hist eta and all phi
      eta_hist[96]->Fill(e);
      eta_hist[ieta]->Fill(e);

      // fill the 2d histo eta, energy and all phi
      energy_eta_hist->Fill(e, ieta);

      // fill 3d histo e_eta_phid
      e_eta_phi->Fill(e, ieta, iphi);
    }

    else if (calotype == LiteCaloEval::HCALOUT)
    {
      // fill the hist with energy data
      // std::cout << ieta << " " <<  iphi  << std::endl;

      if (mode)
      {
        int ket = ieta / 2;

        int llet = ket % 6;
        e *= 0.945 + llet * 0.02;
        int pket = iphi / 4;
        if (pket % 2 == 0)
        {
          e *= 1.03;
        }
      }

      hcal_out_eta_phi[ieta][iphi]->Fill(e);

      hcalout_eta[24]->Fill(e);
      hcalout_eta[ieta]->Fill(e);

      hcalout_energy_eta->Fill(e, ieta);

      hcalout_e_eta_phi->Fill(e, ieta, iphi);
    }

    else if (calotype == LiteCaloEval::HCALIN)
    {
      // fill the hist with energy data

      if (mode)
      {
        int ket = ieta / 2;

        int llet = ket % 6;
        e *= 0.945 + llet * 0.02;
        int pket = iphi / 4;
        if (pket % 2 == 0)
        {
          e *= 1.03;
        }
      }

      hcal_in_eta_phi[ieta][iphi]->Fill(e);

      hcalin_eta[24]->Fill(e);
      hcalin_eta[ieta]->Fill(e);

      hcalin_energy_eta->Fill(e, ieta);

      hcalin_e_eta_phi->Fill(e, ieta, iphi);
    }

    if (!m_UseTowerInfo)
    {
      ++rtiter;
    }

  }  // end of for loop

  _ievent++;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::End(PHCompositeNode * /*topNode*/)
{
  cal_output->cd();
  std::cout << " writing lite calo file" << std::endl;
  cal_output->Write();
  // cout <<" wrote lite calo file" << endl;

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

/// infile histos, outfile is output file name
void LiteCaloEval::Get_Histos(const std::string &infile, const std::string &outfile)
{
  // std::cout << "*** IN GET_HISTOS ***" << std::endl;

  //  std::string outF = outfile;
  //  std::string inF = infile;

  if (infile.empty())
  {
    std::cout << "need infile to run LiteCaloEval::Get_Histos" << std::endl;
    exit(0);
  }

  if (!outfile.empty())
  {
    std::string ts = "cp ";
    ts += infile;
    ts += " ";
    ts += outfile;
    gSystem->Exec(ts.c_str());
    f_temp = new TFile(outfile.c_str(), "UPDATE");  // load the file from the fun4all 1st run
  }

  else
  {
    f_temp = new TFile(infile.c_str());
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

  /// start of eta loop
  for (int i = 0; i < max_ieta + 1; i++)
  {
    std::string b = "eta_" + std::to_string(i);

    if (calotype == LiteCaloEval::HCALOUT)
    {
      b.insert(0, "hcalout_");
    }
    else if (calotype == LiteCaloEval::HCALIN)
    {
      b.insert(0, "hcalin_");
    }

    /// holds the eta slice of histos
    TH1F *heta_temp = (TH1F *) f_temp->Get(b.c_str());

    if (!heta_temp && i == 0)
    {
      std::cout << " warning hist " << b << " not found" << std::endl;
    }

    /// assign the eta slice histo to an array (these arrays are private members in LCE.h)
    eta_hist[i] = heta_temp;

    if (calotype == LiteCaloEval::HCALOUT)
    {
      hcalout_eta[i] = heta_temp;
    }
    else if (calotype == LiteCaloEval::HCALIN)
    {
      hcalin_eta[i] = heta_temp;
    }

    if (!(i < max_ieta))
    {
      continue;
    }

    /// start of phi loop
    for (int j = 0; j < max_iphi; j++)
    {
      /// create string to hold name of individual tower
      std::string hist_name_p = "emc_ieta" + std::to_string(i) + "_phi" + std::to_string(j);

      if (calotype == LiteCaloEval::HCALOUT)
      {
        hist_name_p = "hcal_out_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);
      }
      else if (calotype == LiteCaloEval::HCALIN)
      {
        hist_name_p = "hcal_in_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);
      }

      if (j < 2)
      {
        std::cout << "getting " << hist_name_p << std::endl;
      }

      /// heta_tempp holds tower histogram
      TH1F *heta_tempp = (TH1F *) f_temp->Get(hist_name_p.c_str());

      if (!heta_tempp && i == 0)
      {
        std::cout << " warning hist " << hist_name_p.c_str() << " not found" << std::endl;
      }

      /// assign heta_tempp to array of tower histos
      cemc_hist_eta_phi[i][j] = heta_tempp;

      if (calotype == LiteCaloEval::HCALOUT)
      {
        hcal_out_eta_phi[i][j] = heta_tempp;
      }
      else if (calotype == LiteCaloEval::HCALIN)
      {
        hcal_in_eta_phi[i][j] = heta_tempp;
      }
    }
  }

  // std::cout<< " *** LEAVING GET HISTOS *** "<< std::endl;
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

  float par_value[96] = {0};
  float par_err[96] = {0};
  float eta_value[96] = {0};
  float eta_err[96] = {0};
  //  float rel_err[96];

  if (f_temp)
  {
    f_temp->cd();
  }

  /// fitting - calls the LCE function
  TF1 *f1 = new TF1("myexpo", LCE_fitf, 0.1, 10, 2);
  f1->SetParameters(1.0, 1.0);

  /// if mFS is 1, 10, 20, etc, will make onlyEta = false and will run phi loop below
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
  bool flag_fit_rings = false;

  /// e.x. if mFS = 010 , 010%1000 = 10, 10/100 = 0  -> false, so dont run tower to rings
  /// e.x. if mFS = 110, 110%1000 = 110, 110/100 = 1 -> true, so run to fit eta slice histos fittings
  if (modeFitShifts % 1000 / 100 == 1)
  {
    flag_fit_rings = true;
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

  /// histo to hold the returned fit parameter value
  TH2F *corrPat = new TH2F("corrPat", "", max_ieta, 0, max_ieta, max_iphi, 0, max_iphi);

  if (flag_fit_rings == true)
  {
    if (calotype == LiteCaloEval::CEMC)
    {
      corrPat->SetMinimum(0.95);
      corrPat->SetMaximum(1.03);
    }
    else
    {
      corrPat->SetMinimum(0.92);
      corrPat->SetMaximum(1.05);
    }
  }
  else
  {
    if (calotype == LiteCaloEval::CEMC)
    {
      corrPat->SetMinimum(0.87);
      corrPat->SetMaximum(1.09);
    }
    else
    {
      corrPat->SetMinimum(0.9450);
      corrPat->SetMaximum(1.0764);
    }
  }

  int minbin = 0;
  int maxbin = max_ieta;

  if (m_myminbin > -1)
  {
    minbin = m_myminbin;
  }
  if (m_mymaxbin > -1)
  {
    maxbin = m_mymaxbin;
  }

  /// assign hnewf the eta slice histos.
  TH1F *hnewf = nullptr;

  /// Start of loop for eta only material
  for (int i = minbin; i < maxbin; i++)
  {
    /// a1f, b1f hold the parameter value, parameter error, respectively
    double a1f;
    double b1f;

    /// names for tgraphs, for eta slice histo and graph version of eta slice histo respectively
    std::string myClnm = "newhc_eta" + std::to_string(i);

    std::string myClnm2 = "gnewhc_eta" + std::to_string(i);

    /// dummy indexing
    int iik = i;

    // /// assign hnewf the eta slice histos.
    // TH1F *hnewf = nullptr;
    hnewf = nullptr;

    /// will hold eta slice histo but with certain towers removed, depending on calo type
    TH1F *cleanEtaRef = nullptr;
    std::string cleanEta = "cleanEtaRef_";

    if (calotype == LiteCaloEval::CEMC)
    {
      /// remove TPC
      if (flag_fit_rings == true)
      {
        cleanEta += std::to_string(iik);
        cleanEtaRef = (TH1F *) eta_hist[iik]->Clone(cleanEta.c_str());

        // 4/1/23
        /*
        if( (i >= 2 && i<=93) )
          {
            for(int b = 0; b< 256; b++)
              {
                if( (b >=32 && b <= 34) || (b >= 229 && b <= 231) )
                  {
                    cleanEtaRef->Add((TH1F *)cemc_hist_eta_phi[i][b],-1.0);
                  }
              }
              }*/

        // rebin for run62 to match v9 binning
        cleanEtaRef->Rebin(5);

      }  /// end flag

      else
      {
        hnewf = (TH1F *) ref_lce->eta_hist[iik]->Clone(myClnm.c_str());
        hnewf->Rebin(5);
      }
    }

    else if (calotype == LiteCaloEval::HCALOUT)
    {
      if (flag_fit_rings == true)
      {
        /// remove chimney only from eta bins 0-4, phi bins 14-19
        cleanEta += std::to_string(iik);
        cleanEtaRef = (TH1F *) hcalout_eta[iik]->Clone(cleanEta.c_str());

        if (i < 4)
        {
          for (int b = 14; b <= 19; b++)
          {
            cleanEtaRef->Add((TH1F *) hcal_out_eta_phi[i][b], -1.0);
          }
        }
        // rebin for run62
        cleanEtaRef->Rebin(10);
        cleanEtaRef->Rebin(2);

      }  // end flag

      else
      {
        hnewf = (TH1F *) ref_lce->hcalout_eta[iik]->Clone(myClnm.c_str());
        hnewf->Rebin(10);
        hnewf->Rebin(2);
      }
    }

    else if (calotype == LiteCaloEval::HCALIN)
    {
      if (flag_fit_rings == true)
      {
        cleanEta += std::to_string(iik);
        cleanEtaRef = (TH1F *) hcalin_eta[iik]->Clone(cleanEta.c_str());

        // rebin for run62
        cleanEtaRef->Rebin(10);
        cleanEtaRef->Rebin(10);
      }
      else
      {
        hnewf = (TH1F *) ref_lce->hcalin_eta[iik]->Clone(myClnm.c_str());
        hnewf->Rebin(10);
        hnewf->Rebin(10);
      }
    }

    if (flag_fit_rings == true)
    {
      cleanEtaRef->Smooth(nsmooth);
      LCE_grff = new TGraph(cleanEtaRef);
    }
    else
    {
      /// create graph from the eta slice histo that is sent to fit function
      hnewf->Smooth(nsmooth);
      LCE_grff = new TGraph(hnewf);
    }

    /// this function will be used to fit eta slice histos below
    TF1 *f2f = nullptr;

    /// Fit the eta slices w/ our user defined f'n and then get the tf1
    if (calotype == LiteCaloEval::CEMC)
    {
      eta_hist[i]->Rebin(5);

      if (nsmooth > 1)
      {
        eta_hist[i]->Smooth(nsmooth);
      }

      std::cout << "fitting " << i << std::endl;
      eta_hist[i]->Fit("myexpo", "L", "", fitmin, fitmax);

      f2f = (TF1 *) eta_hist[i]->GetFunction("myexpo");
    }

    else if (calotype == LiteCaloEval::HCALOUT)
    {
      // added 3/3/23 to get the same bin width as previous bin size. Was tested on command line to achieve same bin width as before.
      // old bin width was 0.01, new bin width is 0.0005. So rebin(10) and then rebin(2) worked to get 0.01
      hcalout_eta[i]->Rebin(10);
      hcalout_eta[i]->Rebin(2);

      if (nsmooth > 1)
      {
        hcalout_eta[i]->Smooth(nsmooth);
      }

      hcalout_eta[i]->Fit("myexpo", "L", "", fitmin, fitmax);

      f2f = (TF1 *) hcalout_eta[i]->GetFunction("myexpo");
    }

    else if (calotype == LiteCaloEval::HCALIN)
    {
      // rebinning for run62 only
      hcalin_eta[i]->Rebin(10);
      hcalin_eta[i]->Rebin(10);

      if (nsmooth > 1)
      {
        hcalin_eta[i]->Smooth(nsmooth);
      }

      hcalin_eta[i]->Fit("myexpo", "L", "", fitmin, fitmax);

      f2f = (TF1 *) hcalin_eta[i]->GetFunction("myexpo");
    }

    /// This graph (name is gnewhc_eta) will have points by evaluating
    /// the function (f2f->Eval) retrieved from the fitting above
    TGraph *grT = new TGraph(1000);
    grT->SetName(myClnm2.c_str());

    for (int jjk = 0; jjk < 1000; jjk++)
    {
      float xjj = fitmin + jjk * (fitmax - fitmin) / 1000.0;
      grT->SetPoint(jjk, xjj, f2f->Eval(xjj));
    }

    grT->Write();

    /// Point of this block is to make Fit1_etaout graph

    /// get 'first' parameter, or really the second element bc of indexing
    a1f = f2f->GetParameter(1);
    b1f = f2f->GetParError(1);

    // assign retreived parameter values for graphing
    par_value[i] = a1f;
    par_err[i] = b1f;
    eta_value[i] = i;
    eta_err[i] = 0.01;
    // rel_err[i] = par_err[i] / par_value[i];

    /// This if statement enables the phi (tower) material to execute
    if (onlyEta == true)
    {
      continue;
    }

    /****************************************************
      This is the nested forloop to start tower material
     ****************************************************/

    for (int j = 0; j < max_iphi; j++)
    {
      // 4/1/23
      /*
      /// skip over iterations associated with TPC regions / chimney
      if(flag_fit_rings == true)
        {

          if(calotype == LiteCaloEval::CEMC)
            {

              if(i >= 2 && i <= 93)
                {
                  if( (j >=32 && j <= 34) || (j >=229 && j <= 231) )
                    {
                      continue;
                    }
                }

            }
          else if(calotype == LiteCaloEval::HCALOUT)
            {
              /// skip over chimney
              if( (i >=0 && i <=3) && (j>=14 && j <= 19))
                {
                  continue;
                }
            }

        }//end flag
          */

      /// names of tower histo and graph version of same histo
      std::string myClnmp = "newhc_eta" + std::to_string(1000 * (i + 2) + j);

      std::string myClnm2p = "gnewhc_eta" + std::to_string(1000 * (i + 2) + j);

      if (j == 0)
      {
        std::cout << " making fitting  " << myClnmp << std::endl;
      }

      /// histo to hold tower from non modified root files
      TH1F *hnewfp = nullptr;

      if (calotype == LiteCaloEval::CEMC)
      {
        cemc_hist_eta_phi[i][j]->Rebin(5);

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) cemc_hist_eta_phi[i][j], -1.0);
        }

        else
        {
          hnewfp = (TH1F *) ref_lce->cemc_hist_eta_phi[i][j]->Clone(myClnmp.c_str());
          hnewfp->Rebin(5);
        }
      }

      else if (calotype == LiteCaloEval::HCALOUT)
      {
        // for run62 added 3/3/23 to get the same bin width as previous bins. Was tested on command line to achieve same bin width as before.
        // old bin width was 0.01, new bin width is 0.0005. So rebin(10) and then rebin(2) worked to get 0.01
        hcal_out_eta_phi[i][j]->Rebin(10);
        hcal_out_eta_phi[i][j]->Rebin(2);

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) hcal_out_eta_phi[i][j], -1.0);
        }

        else
        {
          hnewfp = (TH1F *) ref_lce->hcal_out_eta_phi[i][j]->Clone(myClnmp.c_str());
          hnewfp->Rebin(10);
          hnewfp->Rebin(2);
        }
      }

      else if (calotype == LiteCaloEval::HCALIN)
      {
        // rebin tower for run62
        hcal_in_eta_phi[i][j]->Rebin(10);
        hcal_in_eta_phi[i][j]->Rebin(10);

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) hcal_in_eta_phi[i][j], -1.0);
        }
        else
        {
          hnewfp = (TH1F *) ref_lce->hcal_in_eta_phi[i][j]->Clone(myClnmp.c_str());
          hnewfp->Rebin(10);
          hnewfp->Rebin(10);
        }
      }

      if (j < 2)
      {
        std::cout << "got neweff phi eta  ... fitting w/ fit min,max: " << fitmin
                  << " " << fitmax << std::endl;
      }

      /*
       If false make tgraph out of tower and send to be fit with fit fuction.
       If true the towers are actually fit against the eta ring. This happends bc "myexpo" uses a tgraph
       to fit, and if we dont make LCE_grff for a tower, myexpo uses the latest version of LCE_grff, which is an eta slice tgraph
      */
      if (flag_fit_rings == false)
      {
        hnewfp->Smooth(nsmooth);
        LCE_grff = new TGraph(hnewfp);
      }

      if (flag_fit_rings == true)
      {
        // cleanEtaRef->Smooth(nsmooth);
        LCE_grff = new TGraph(cleanEtaRef);
      }

      /// make tf1 that will hold the resulting fit from myexpo on towers
      TF1 *f2f2 = nullptr;

      if (calotype == LiteCaloEval::CEMC)
      {
        double scaleP0 = cemc_hist_eta_phi[i][j]->Integral(cemc_hist_eta_phi[i][j]->FindBin(fitmin), cemc_hist_eta_phi[i][j]->FindBin(fitmax));

        if (flag_fit_rings == 1)
        {
          scaleP0 /= cleanEtaRef->Integral(cleanEtaRef->FindBin(fitmin), cleanEtaRef->FindBin(fitmax));
        }
        else
        {
          scaleP0 /= hnewfp->Integral(hnewfp->FindBin(fitmin), hnewfp->FindBin(fitmax));
        }

        f1->SetParameters(scaleP0, 1.0);

        if (j < 2)
        {
          cemc_hist_eta_phi[i][j]->Fit("myexpo", "L", "", fitmin, fitmax);
        }

        else
        {
          cemc_hist_eta_phi[i][j]->Fit("myexpo", "LQ", "", fitmin, fitmax);
        }

        f2f2 = (TF1 *) cemc_hist_eta_phi[i][j]->GetFunction("myexpo");

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) cemc_hist_eta_phi[i][j], 1.0);
        }
      }

      else if (calotype == LiteCaloEval::HCALOUT)
      {
        double scaleP0 = hcal_out_eta_phi[i][j]->Integral(hcal_out_eta_phi[i][j]->FindBin(fitmin), hcal_out_eta_phi[i][j]->FindBin(fitmax));

        if (flag_fit_rings == 1)
        {
          scaleP0 /= cleanEtaRef->Integral(cleanEtaRef->FindBin(fitmin), cleanEtaRef->FindBin(fitmax));
        }
        else
        {
          scaleP0 /= hnewfp->Integral(hnewfp->FindBin(fitmin), hnewfp->FindBin(fitmax));
        }

        /// reset the "myexpo" parameters
        f1->SetParameters(scaleP0, 1.0);

        if (j < 2)
        {
          hcal_out_eta_phi[i][j]->Fit("myexpo", "L", "", fitmin, fitmax);
        }

        else
        {
          hcal_out_eta_phi[i][j]->Fit("myexpo", "LQ", "", fitmin, fitmax);
        }

        f2f2 = (TF1 *) hcal_out_eta_phi[i][j]->GetFunction("myexpo");

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) hcal_out_eta_phi[i][j], 1.0);
        }
      }

      else if (calotype == LiteCaloEval::HCALIN)
      {
        double scaleP0 = hcal_in_eta_phi[i][j]->Integral(hcal_in_eta_phi[i][j]->FindBin(fitmin), hcal_in_eta_phi[i][j]->FindBin(fitmax));

        if (flag_fit_rings == 1)
        {
          scaleP0 /= cleanEtaRef->Integral(cleanEtaRef->FindBin(fitmin), cleanEtaRef->FindBin(fitmax));
        }

        else
        {
          scaleP0 /= hnewfp->Integral(hnewfp->FindBin(fitmin), hnewfp->FindBin(fitmax));
        }

        f1->SetParameters(scaleP0, 1.0);

        if (j < 2)
        {
          hcal_in_eta_phi[i][j]->Fit("myexpo", "L", "", fitmin, fitmax);
        }
        else
        {
          hcal_in_eta_phi[i][j]->Fit("myexpo", "LQ", "", fitmin, fitmax);
        }

        f2f2 = (TF1 *) hcal_in_eta_phi[i][j]->GetFunction("myexpo");

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) hcal_in_eta_phi[i][j], 1.0);
        }
      }

      /// make graph and set points from evaluating f2f2
      TGraph *local_grT = new TGraph(1000);
      local_grT->SetName(myClnm2p.c_str());

      for (int jjk = 0; jjk < 1000; jjk++)
      {
        float xjj = fitmin + jjk * (fitmax - fitmin) / 1000.0;
        // grT->SetPoint(jjk,xjj,f2f->Eval(xjj));
        local_grT->SetPoint(jjk, xjj, f2f2->Eval(xjj));
      }

      local_grT->Write();

      float a1fp = f2f2->GetParameter(1);
      float b1fp = f2f2->GetParError(1);
      // float b1fp = f2f->GetParError(1);

      /// fill corrpat, which has returned fit values from towers
      /// e.x. if you want to view the histo that is, say, eta=12, phi = 1,
      /// you really need to draw the histo with eta=11,phi=0 bc of the i+1, j+1
      corrPat->SetBinContent(i + 1, j + 1, 1 / a1fp);
      corrPat->SetBinError(i + 1, j + 1, 1 / b1fp);

    }  // end of inner forloop (phi)

  }  // end of outter forloop (eta)

  // create graph that plots eta slice par values
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
  }

  f_temp->Close();

  // std::cout << "LEAVING REL SHIFTS" << std::endl;
}
