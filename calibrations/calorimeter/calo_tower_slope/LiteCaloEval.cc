#include "LiteCaloEval.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <ffarawobjects/Gl1Packet.h>
#include <calotrigger/TriggerAnalyzer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>


#include <RtypesCore.h>  // for Double_t
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <boost/format.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

class RawTowerGeom;

/// This function is used for the histo fitting process. x is a 1d array that holds xaxis values.
/// par is an array of 1d array of parameters we set our fit function to. So p[0] = p[1] = 1 unless otherwise noted
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
  std::cout << "In InitRun " << std::endl;

  // just quit if we forgot to set the calorimeter type
  if (calotype == LiteCaloEval::NONE)
  {
    std::cout << Name() << ": No Calo Type set" << std::endl;
    gSystem->Exit(1);
  }

  _ievent = 0;

  trigAna = new TriggerAnalyzer();

  cal_output = new TFile(_filename.c_str(), "RECREATE");

  h_event = new TH1F("h_event","",1,0,1);

  if (calotype == LiteCaloEval::HCALIN)
  {
    hcalin_energy_eta = new TH2F("hcalin_energy_eta", "hcalin energy eta", 100, 0, 10, 24, -0.5, 23.5);
    hcalin_e_eta_phi = new TH3F("hcalin_e_eta_phi", "hcalin e eta phi", 60, 0, 6, 24, -0.5, 23.5, 64, -0.5, 63.5);

    /// create tower histos
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_in_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_in_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_in_energy", 40000, 0, 4);
        hcal_in_eta_phi[i][j]->SetXTitle("Energy [GeV]");
      }
    }

    // create eta slice histos
    for (int i = 0; i < 25; i++)
    {
      std::string hist_name = "hcalin_eta_" + std::to_string(i);

      if (i < 24)
      {
        hcalin_eta[i] = new TH1F(hist_name.c_str(), "hcalin eta's", 4000, 0, 4.);
        hcalin_eta[i]->SetXTitle("Energy [GeV]");
      }
      else
      {
        hcalin_eta[i] = new TH1F(hist_name.c_str(), "hcalin eta's", 40000, 0, 4.);
        hcalin_eta[i]->SetXTitle("Energy [GeV]");
      }
    }
  }

  else if (calotype == LiteCaloEval::HCALOUT)
  {
    hcalout_energy_eta = new TH2F("hcalout_energy_eta", "hcalout energy eta", 100, 0, 10, 24, 0.5, 23.5);
    hcalout_e_eta_phi = new TH3F("hcalout_e_eta_phi", "hcalout e eta phi", 100, 0, 10, 24, -0.5, 23.5, 64, -0.5, 63.5);

    /// create tower histos
    for (int i = 0; i < 24; i++)
    {
      for (int j = 0; j < 64; j++)
      {
        std::string hist_name = "hcal_out_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);

        hcal_out_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hcal_out energy", 10000, 0, 10);
        hcal_out_eta_phi[i][j]->SetXTitle("Energy [GeV]");
      }
    }

    /// create eta slice histos
    for (int i = 0; i < 25; i++)
    {
      std::string hist_name = "hcalout_eta_" + std::to_string(i);
      if (i < 24)
      {
        hcalout_eta[i] = new TH1F(hist_name.c_str(), "hcalout eta's", 10000, 0, 10);
        hcalout_eta[i]->SetXTitle("Energy [GeV]");
      }
      else
      {
        hcalout_eta[i] = new TH1F(hist_name.c_str(), "hcalout eta's", 100000, 0, 10);
        hcalout_eta[i]->SetXTitle("Energy [GeV]");
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

        cemc_hist_eta_phi[i][j] = new TH1F(hist_name.c_str(), "Hist_ieta_phi_leaf(e)", 400, 0, 2);
        cemc_hist_eta_phi[i][j]->SetXTitle("Energy [GeV]");
      }
    }

    // create eta slice histos
    for (int i = 0; i < 97; i++)
    {
      gStyle->SetOptFit(1);
      std::string b = "eta_" + std::to_string(i);
      eta_hist[i] = new TH1F(b.c_str(), "eta and all phi's", 400, 0, 2);
      eta_hist[i]->SetXTitle("Energy [GeV]");
    }

    // make 2d histo
    energy_eta_hist = new TH2F("energy_eta_hist", "energy eta and all phi", 100, 0, 10, 96, -0.5, 95.5);

    // make 3d histo
    e_eta_phi = new TH3F("e_eta_phi", "e v eta v phi", 100, 0, 10, 96, -0.5, 95.5, 256, -0.5, 255.5);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LiteCaloEval::process_event(PHCompositeNode *topNode)
{
  if (_ievent % 100 == 0)
  {
    std::cout << "LiteCaloEval::process_event(PHCompositeNode *topNode) Processing Event " << _ievent << std::endl;
  }

  //--------------------------- trigger and GL1-------------------------------//
  trigAna->decodeTriggers(topNode);

  if (reqMinBias && trigAna->didTriggerFire(12) == false)
  {
    _ievent++;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  h_event->Fill(0);

  //---------------------------- Get geometry -------------------------------------//
  // raw tower container
  std::string towernode = "TOWER_CALIB_" + _caloname;
  RawTowerContainer *towers = nullptr;
  RawTowerGeomContainer *towergeom = nullptr;

  // for tower energy
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

  if (mode && _ievent < 2)
  {
    std::cout << "mode is set " << std::endl;
  }

  TowerInfoContainer *towerinfos = nullptr;

  // if using towerinfo create a tower object
  if (m_UseTowerInfo)
  {
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

  for (unsigned int channel = 0; channel < nchannels; channel++)
  {
    if (!m_UseTowerInfo && rtiter == begin_end.second)
    {
      std::cout << "ERROR: Recheck iteration process with rtiter variable" << std::endl;
    }

    float e = -1.0;  // tower energy
    unsigned int ieta = 99999;
    unsigned int iphi = 99999;

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
      // mode variable is used for simulations only. Creates a decalibration in energy
      if (mode)
      {
        int ket = ieta / 8;
        int llet = ket % 6;
        int pket = iphi / 8;
        int ppkket = pket % 3;

        e *= 0.88 + llet * 0.04 - 0.01 + 0.01 * ppkket;
      }

      cemc_hist_eta_phi[ieta][iphi]->Fill(e);

      eta_hist[96]->Fill(e);

      eta_hist[ieta]->Fill(e);

      energy_eta_hist->Fill(e, ieta);

      e_eta_phi->Fill(e, ieta, iphi);
    }

    else if (calotype == LiteCaloEval::HCALOUT)
    {
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

  return Fun4AllReturnCodes::EVENT_OK;
}

/// infile histos, outfile is output file name
void LiteCaloEval::Get_Histos(const std::string &infile, const std::string &outfile)
{
  std::cout << "Getting histograms... " << std::endl;

  if (infile.empty())
  {
    std::cout << "Need an input file to run LiteCaloEval::Get_Histos" << std::endl;
    exit(0);
  }

  if (!outfile.empty())
  {
    std::string ts = "cp ";
    ts += infile;
    ts += " ";
    ts += outfile;
    gSystem->Exec(ts.c_str());
    f_temp = new TFile(outfile.c_str(), "UPDATE");
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

    if (i == 0)
    {
      std::cout << "Old bin width for eta slice " << heta_temp->GetBinWidth(1) << std::endl;
    }

    // rebin histogram
    heta_temp->Rebin(binwidth / heta_temp->GetBinWidth(1));

    if (i == 0)
    {
      std::cout << "New bin width for eta slice " << heta_temp->GetBinWidth(1) << std::endl;
    }

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

      /// heta_tempp holds tower histogram
      TH1F *heta_tempp = (TH1F *) f_temp->Get(hist_name_p.c_str());

      if (i == 0 && j == 0)
      {
        std::cout << "Old bin width for towers " << heta_tempp->GetBinWidth(1) << std::endl;
      }

      // rebin histogram
      heta_tempp->Rebin(binwidth / heta_tempp->GetBinWidth(1));

      if (i == 0 && j == 0)
      {
        std::cout << "New bin width for towers " << heta_tempp->GetBinWidth(1) << std::endl;
      }

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
    }//phi loop
  }//eta loop

  /*
  f_temp->Close();
  f_temp = nullptr;
  delete f_temp;
  */

  std::cout << "Grabbed all histograms." << std::endl;

}//end Get_Histos f'n

void LiteCaloEval::FitRelativeShifts(LiteCaloEval *ref_lce, int modeFitShifts)
{
  bool onlyEta = false;  // will determine if you only run over eta slices or not

  if (fitmin < 0.001)
  {
    fitmin = 0.15;
  }
  if (fitmax < 0.001)
  {
    fitmax = 1.3;
  }

  float par_value[96] = {0};  // gain value used only in the outer for loop
  float par_err[96] = {0};    // error on the gain
  float eta_value[96] = {0};  // ieta
  float eta_err[96] = {0};    // ieta err. just will be zero.

  if (f_temp)
  {
    f_temp->cd();
  }

  /// fitting - calls the LCE_fitf function at beginning of module
  TF1 *f1 = new TF1("myexpo", LCE_fitf, 0.1, 10, 2);
  f1->SetParameters(1.0, 1.0);

  /// looks at 1's palce of mFS parameter, if onlyEta = false will run phi loop (individual towers) below
  if (modeFitShifts % 10 == 1)
  {
    onlyEta = true;
    std::cout << "Running only eta slice mode....." << std::endl;
  }

  /// nsmooth is an automatic smoothing of histos
  int nsmooth = 1;

  /// if want more smoothing, look at tens place in mFS. If mFS=10, then addSmooth=1, and then there will be 2 total smoothings
  int addSmooth = modeFitShifts % 100 / 10;
  if (addSmooth)
  {
    nsmooth += addSmooth;
  }

  /// flag that will run fits at tower level to get shifts, or will run at ring level and fit towers to eta slice histo
  bool flag_fit_rings = false;

  // look at 100's place of mFS. Runs the eta slice flattening mode. If false, runs gain tracing mode
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
  corrPat->SetXTitle("#eta bin");
  corrPat->SetYTitle("#phi bin");

  TH2F *h2_failQA = new TH2F("h2_failQA", "", max_ieta, 0, max_ieta, max_iphi, 0, max_iphi);

  // 1d histo for gain shift values
  TH1F *gainvals = new TH1F("gainvals", "Towerslope Correction Values", 10000, 0, 10);
  gainvals->SetXTitle("Gain Shift Value");
  gainvals->SetYTitle("Counts");

  // 1d histo for gain shift error
  TH1F *h_gainErr = new TH1F("h_gainErr", "Towerslope Corrections Errors", 1000, 0, 1);
  h_gainErr->SetXTitle("error");
  h_gainErr->SetYTitle("Counts");

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

  /// assign hnewf the eta slice histos when running in Gain Trace mode
  TH1F *hnewf = nullptr;

  /// Start of loop for eta slices
  for (int i = minbin; i < maxbin; i++)
  {
    std::cout << "===============================" << std::endl;
    std::cout << "Fitting towers in eta slice " << i << std::endl;

    /// hold the gain value, gain error, respectively from the fit
    double ieta_gain;
    double ieta_gain_err;

    /// names for eta slice histo
    std::string myClnm = "newhc_eta" + std::to_string(i);

    /// dummy indexing for string name
    int iik = i;

    /// will hold eta slice histo clone but with current tower being fitted removed from the eta slice
    TH1F *cleanEtaRef = nullptr;
    std::string cleanEta = "cleanEtaRef_";

    if (calotype == LiteCaloEval::CEMC)
    {
      if (flag_fit_rings == true)
      {
        cleanEta += std::to_string(iik);

        cleanEtaRef = (TH1F *) eta_hist[iik]->Clone(cleanEta.c_str());
      }

      else
      {
        hnewf = (TH1F *) ref_lce->eta_hist[iik]->Clone(myClnm.c_str());
      }
    }

    else if (calotype == LiteCaloEval::HCALOUT)
    {
      if (flag_fit_rings == true)
      {
        cleanEta += std::to_string(iik);

        cleanEtaRef = (TH1F *) hcalout_eta[iik]->Clone(cleanEta.c_str());

        // remove towers from eta slice reference associated with the chimney (i < 4) and support ring in high eta region(i>19)
        if (i < 4 || i > 19)
        {
          for (int phiCH = 14; phiCH < 20; phiCH++)
          {
            cleanEtaRef->Add((TH1F *) hcal_out_eta_phi[i][phiCH], -1.0);
          }
        }
      }

      else
      {
        hnewf = (TH1F *) ref_lce->hcalout_eta[iik]->Clone(myClnm.c_str());
      }
    }

    else if (calotype == LiteCaloEval::HCALIN)
    {
      if (flag_fit_rings == true)
      {
        cleanEta += std::to_string(iik);

        cleanEtaRef = (TH1F *) hcalin_eta[iik]->Clone(cleanEta.c_str());
      }

      else
      {
        hnewf = (TH1F *) ref_lce->hcalin_eta[iik]->Clone(myClnm.c_str());
      }
    }

    if (flag_fit_rings == true)
    {
      cleanEtaRef->Smooth(nsmooth);

      LCE_grff = new TGraph(cleanEtaRef);
    }

    else
    {
      hnewf->Smooth(nsmooth);

      LCE_grff = new TGraph(hnewf);
    }

    /// this function will be used to fit eta slice histos
    TF1 *f2f = nullptr;

    /// Fit the eta slices w/ our user defined f'n and then get the tf1
    if (calotype == LiteCaloEval::CEMC)
    {
      if (nsmooth > 1)
      {
        eta_hist[i]->Smooth(nsmooth);
      }

      eta_hist[i]->Fit("myexpo", "Q", "", fitmin, fitmax);

      f2f = (TF1 *) eta_hist[i]->GetFunction("myexpo");
    }

    else if (calotype == LiteCaloEval::HCALOUT)
    {
      if (nsmooth > 1)
      {
        hcalout_eta[i]->Smooth(nsmooth);
      }

      hcalout_eta[i]->Fit("myexpo", "Q", "", fitmin, fitmax);

      f2f = (TF1 *) hcalout_eta[i]->GetFunction("myexpo");
    }

    else if (calotype == LiteCaloEval::HCALIN)
    {
      if (nsmooth > 1)
      {
        hcalin_eta[i]->Smooth(nsmooth);
      }

      hcalin_eta[i]->Fit("myexpo", "Q", "", fitmin, fitmax);

      f2f = (TF1 *) hcalin_eta[i]->GetFunction("myexpo");
    }

    if(!f2f)
      {
	std::cout << "Warning, f2f is null!" << std::endl;
	exit(-1);
      }
    else
      {
	ieta_gain = f2f->GetParameter(1);
      }

    ieta_gain_err = f2f->GetParError(1);

    par_value[i] = ieta_gain;
    par_err[i] = ieta_gain_err;
    eta_value[i] = i;
    eta_err[i] = 0.0;

    /// This if statement enables the tower fits to run
    if (onlyEta == true)
    {
      continue;
    }

    /////////////////////
    // QA -- Blair
    if (doQA)
    {
      if (flag_fit_rings == true)
      {
        if (calotype == LiteCaloEval::CEMC)
        {
          for (int j = 0; j < max_iphi; j++)
          {
            cleanEtaRef->Add((TH1F *) cemc_hist_eta_phi[i][j], -1.0);

            bool qa_res = spec_QA(cemc_hist_eta_phi[i][j], cleanEtaRef, true);
            if (qa_res == false)
            {
              std::cout << "failed QA " << i << "," << j << std::endl;
              cemc_hist_eta_phi[i][j]->Reset();
              h2_failQA->Fill(i, j);
            }
            else
            {
              cleanEtaRef->Add((TH1F *) cemc_hist_eta_phi[i][j], 1.0);
            }
          }
        }
      }
    }
    /****************************************************
    This is the nested forloop to start tower material
    ****************************************************/

    for (int j = 0; j < max_iphi; j++)
    {
      /// names of tower histo for cloning. used in gain trace mode
      std::string myClnmp = "newhc_eta" + std::to_string(1000 * (i + 2) + j);

      /// histo to hold tower clone
      TH1F *hnewfp = nullptr;

      // check to see if tower in question is part of the chimney/high eta support ring

      bool _isChimney = chk_isChimney(i, j);

      if (calotype == LiteCaloEval::CEMC)
      {
        // skip tower if there are no entries
        if (!(cemc_hist_eta_phi[i][j]->GetEntries()))
        {
          std::cout << "No entries in EMCAL tower histogram (" << i << "," << j << "). Skipping fitting." << std::endl;
          continue;
        }

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) cemc_hist_eta_phi[i][j], -1.0);
        }
        else
        {
          hnewfp = (TH1F *) ref_lce->cemc_hist_eta_phi[i][j]->Clone(myClnmp.c_str());
        }
      }

      else if (calotype == LiteCaloEval::HCALOUT)
      {
        // skip tower if there are no entries
        if (!(hcal_out_eta_phi[i][j]->GetEntries()))
        {
          std::cout << "No entries in OHCAL tower histogram (" << i << "," << j << "). Skipping fitting." << std::endl;
          continue;
        }

        if (flag_fit_rings == true)
        {
          // check to see if tower is part of chimney/high eta support. If it isnt, proceed to remove that tower from eta ref. This is to ensure we dont remove towers in these regions  twice from eta ref bc they already were in L717

          if (!_isChimney)
          {
            cleanEtaRef->Add((TH1F *) hcal_out_eta_phi[i][j], -1.0);
          }
        }

        else
        {
          hnewfp = (TH1F *) ref_lce->hcal_out_eta_phi[i][j]->Clone(myClnmp.c_str());
        }
      }

      else if (calotype == LiteCaloEval::HCALIN)
      {
        // skip tower if there are no entries
        if (!(hcal_in_eta_phi[i][j]->GetEntries()))
        {
          std::cout << "No entries in IHCAL tower histogram(" << i << "," << j << "). Skipping fitting." << std::endl;
          continue;
        }

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) hcal_in_eta_phi[i][j], -1.0);
        }
        else
        {
          hnewfp = (TH1F *) ref_lce->hcal_in_eta_phi[i][j]->Clone(myClnmp.c_str());
        }
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

      /// make tf1 that will hold the resulting fit of towers
      TF1 *f2f2 = nullptr;

      // need to scale the reference histogram to allow it to start at an amplitude similar/at tower that is to be fit
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

        cemc_hist_eta_phi[i][j]->Fit("myexpo", "Q", "", fitmin, fitmax);

        f2f2 = (TF1 *) cemc_hist_eta_phi[i][j]->GetFunction("myexpo");

        // add back the just fitted tower to the eta slice reference
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

        f1->SetParameters(scaleP0, 1.0);

        hcal_out_eta_phi[i][j]->Fit("myexpo", "Q", "", fitmin, fitmax);

        f2f2 = (TF1 *) hcal_out_eta_phi[i][j]->GetFunction("myexpo");

        if (flag_fit_rings == true)
        {
          if (!_isChimney)
          {
            cleanEtaRef->Add((TH1F *) hcal_out_eta_phi[i][j], 1.0);
          }
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

        hcal_in_eta_phi[i][j]->Fit("myexpo", "Q", "", fitmin, fitmax);

        f2f2 = (TF1 *) hcal_in_eta_phi[i][j]->GetFunction("myexpo");

        if (flag_fit_rings == true)
        {
          cleanEtaRef->Add((TH1F *) hcal_in_eta_phi[i][j], 1.0);
        }
      }

      float correction = f2f2->GetParameter(1);

      float corrErr = f2f2->GetParError(1);

      double errProp = corrErr / (correction * correction);

      /// fill corrpat, which has returned fit values from towers
      /// e.x. if you want to view the histo that is, say, eta=12, phi = 1,
      /// you really need to draw the histo with eta=11,phi=0 bc of the i+1, j+1

      corrPat->SetBinContent(i + 1, j + 1, 1 / correction);

      corrPat->SetBinError(i + 1, j + 1, errProp);

      gainvals->Fill(1.0 / correction);

      h_gainErr->Fill(errProp);

    }  // end of inner forloop (phi)

  }  // end of outter forloop (eta)

  // create graph that plots eta slice par values (this is only when looping over eta slices - not towers)
  TGraphErrors g1(max_ieta, eta_value, par_value, eta_err, par_err);
  g1.SetTitle("fitted shifts; eta; p1");
  g1.SetMarkerStyle(20);
  g1.Draw("ap");
  g1.SetName("Fit1_etaout");
  g1.Write();

  corrPat->Write();
  h2_failQA->Write();
  gainvals->Write();
  h_gainErr->Write();

 

  if (f_temp)
  {
    f_temp->Write();
  }

  f_temp->Close();
}

bool LiteCaloEval::chk_isChimney(int ieta, int iphi)
{
  if ((ieta < 4 || ieta > 19) && (iphi > 13 && iphi < 20))
  {
    return true;
  }

  return false;
}

float LiteCaloEval::spec_QA(TH1 *h_spec, TH1 *h_ref, bool /*retFloat*/)
{
  float norm = h_spec->Integral(h_spec->FindBin(0.4, 1.1), h_spec->FindBin(1.0));
  if (norm == 0)
  {
    std::cout << "no data for QA" << std::endl;
    return false;
  }
  TF1 *f_exp = new TF1("f_exp", "[0]*TMath::Exp(-1.0*x*[1])", 0.5, 1);
  h_spec->Fit(f_exp, "QNR", "", 0.4, 1.);
  float spec_p1 = f_exp->GetParameter(1);
  h_ref->Fit(f_exp, "QNR", "", 0.4, 1.);
  float ref_p1 = f_exp->GetParameter(1);
  float calib = ref_p1 / spec_p1;
  delete f_exp;
  return calib;
}

bool LiteCaloEval::spec_QA(TH1 *h_spec, TH1 *h_ref)
{
  float calib = spec_QA(h_spec, h_ref, true);
  if (calib > 1.2 || calib < 0.8)
  {
    return false;
  }
  else
  {
    return true;
  }
}

void LiteCaloEval::plot_cemc(const std::string &path)
{
  TH2F *h_fail = new TH2F("h_fail", "", 96, 0, 96, 256, 0, 256);
  TH2F *h_corrPat = new TH2F("corrPat", "", 96, 0, 96, 256, 0, 256);
  TH2F *h_hits = new TH2F("h_hits", "", 96, 0, 96, 256, 0, 256);
  TH2F *h_expCalib = new TH2F("h_expCalib", "", 96, 0, 96, 256, 0, 256);
  h_hits->SetXTitle("#it{#eta}_{i}");
  h_hits->SetYTitle("#it{#phi}_{i}");
  TH2F *h_hits_clean = new TH2F("h_hits_clean", "", 96, 0, 96, 256, 0, 256);
  h_hits_clean->SetXTitle("#it{#eta}_{i}");
  h_hits_clean->SetYTitle("#it{#phi}_{i}");

  for (int ieta = 0; ieta < 96; ieta++)
  {
    for (int iphi = 0; iphi < 256; iphi++)
    {
      h_corrPat->SetBinContent(ieta + 1, iphi + 1, 1);
      cemc_hist_eta_phi[ieta][iphi]->Rebin(5);
      // record hits
      int binhit = cemc_hist_eta_phi[ieta][iphi]->FindBin(0.25);
      int binMax = cemc_hist_eta_phi[ieta][iphi]->GetNbinsX();
      float hits = cemc_hist_eta_phi[ieta][iphi]->Integral(binhit, binMax + 1);
      h_hits->SetBinContent(ieta + 1, iphi + 1, hits);
      h_hits->SetBinError(ieta + 1, iphi + 1, sqrt(hits));
      h_hits_clean->SetBinContent(ieta + 1, iphi + 1, hits);
      h_hits_clean->SetBinError(ieta + 1, iphi + 1, sqrt(hits));
      int bin = cemc_hist_eta_phi[ieta][iphi]->FindBin(0.2);
      if (cemc_hist_eta_phi[ieta][iphi]->GetBinContent(bin) == 0)
      {
        continue;
      }
      float calib = spec_QA(cemc_hist_eta_phi[ieta][iphi], eta_hist[ieta], true);
      h_expCalib->SetBinContent(ieta + 1, iphi + 1, calib);
      if (!spec_QA(cemc_hist_eta_phi[ieta][iphi], eta_hist[ieta], true))
      {
        h_fail->Fill(ieta, iphi);
        h_corrPat->SetBinContent(ieta + 1, iphi + 1, 0);
        h_hits_clean->SetBinContent(ieta + 1, iphi + 1, 0);
      }
    }
  }

  // find hot towers
  TH2F *h_hot = new TH2F("h_hot", "", 96, 0, 96, 256, 0, 256);
  TH1F *h1_hits[96];
  TH1F *h1_hits2[96];
  float max = h_hits_clean->GetBinContent(h_hits_clean->GetMaximumBin());
  float min = h_hits_clean->GetMinimum();
  if (min == 0)
  {
    min = 1;
  }

  TCanvas *c3 = new TCanvas("c3", "c3", 1600, 600);

  for (int ie = 0; ie < 96; ie++)
  {
    std::string name1 = "h1_hits" + std::to_string(ie);
    std::string name2 = "h1_hits2_" + std::to_string(ie);

    h1_hits[ie] = new TH1F(name1.c_str(), "", 200, min, max);
    h1_hits2[ie] = new TH1F(name2.c_str(), "", 200, min, max);
    for (int iphi = 0; iphi < 256; iphi++)
    {
      h1_hits[ie]->Fill(h_hits_clean->GetBinContent(ie + 1, iphi + 1));
    }
    float mean = h1_hits[ie]->GetMean();
    float std = h1_hits[ie]->GetStdDev();
    for (int iphi = 0; iphi < 256; iphi++)
    {
      float val = h_hits_clean->GetBinContent(ie + 1, iphi + 1);
      if (fabs(mean - val) / std > 3)
      {
        continue;
      }
      h1_hits2[ie]->Fill(val);
    }
    mean = h1_hits[ie]->GetMean();
    std = h1_hits[ie]->GetStdDev();
    for (int iphi = 0; iphi < 256; iphi++)
    {
      h_hot->SetBinContent(ie + 1, iphi + 1, 0);
      float val = h_hits_clean->GetBinContent(ie + 1, iphi + 1);
      if (val == 0)
      {
        continue;
      }
      if (fabs(mean - val) / std > 5)
      {
        h_hot->SetBinContent(ie + 1, iphi + 1, 1);
        h_corrPat->SetBinContent(ie + 1, iphi + 1, 0);
        h_hits_clean->SetBinContent(ie + 1, iphi + 1, 0);
      }
    }
    h1_hits[ie]->Write();
    h1_hits2[ie]->Write();
    if (ie < 90)
    {
      continue;
    }
    h1_hits[ie]->Draw("hist same");
    h1_hits[ie]->GetYaxis()->SetRangeUser(0.1, 100);
    h1_hits[ie]->GetXaxis()->SetTitle("Number of hits");
    h1_hits[ie]->GetXaxis()->SetNdivisions(505);
    h1_hits2[ie]->Draw("hist same");
    h1_hits2[ie]->SetLineColor(kBlue);
  }
  c3->SaveAs((path + "/hits.pdf").c_str());

  TLatex latex;
  latex.SetTextFont(42);
  latex.SetTextSize(0.04);
  latex.SetTextAlign(12);  // Align to the left
  latex.SetNDC();

  int ieta = 0;
  for (int iplot = 0; iplot < 10; iplot++)
  {
    TCanvas *c1 = new TCanvas("c1", "c1", 400, 600);

    for (int ie = 0; ie < 10; ie++)
    {
      if (ieta > 95)
      {
        continue;
      }
      for (int iphi = 0; iphi < 256; iphi++)
      {
        // scale hist
        int bin = cemc_hist_eta_phi[ieta][iphi]->FindBin(0.2);
        if (cemc_hist_eta_phi[ieta][iphi]->GetBinContent(bin) == 0)
        {
          continue;
        }
        cemc_hist_eta_phi[ieta][iphi]->Scale(1.0 / cemc_hist_eta_phi[ieta][iphi]->GetBinContent(bin) * pow(10, ie));
        // plot
        cemc_hist_eta_phi[ieta][iphi]->Draw("same hist");
        cemc_hist_eta_phi[ieta][iphi]->SetXTitle("pre-calib tower [GeV]");
        cemc_hist_eta_phi[ieta][iphi]->SetYTitle("N_{twr}/N_{twr}(E=0.2 GeV) x10^{ieta}");
        cemc_hist_eta_phi[ieta][iphi]->GetXaxis()->SetRangeUser(0, 1);
        cemc_hist_eta_phi[ieta][iphi]->GetYaxis()->SetRangeUser(0.01, 1e11);
        if (h_corrPat->GetBinContent(ieta + 1, iphi + 1) == 0)
        {
          cemc_hist_eta_phi[ieta][iphi]->SetLineColor(kBlue);
          h_fail->Fill(ieta, iphi);
          h_hits_clean->SetBinContent(ieta + 1, iphi + 1, 0);
        }
        latex.SetTextSize(0.02);
        latex.SetTextColor(kBlack);
        std::string result = "#eta_{" + std::to_string(ieta) + "}#times10^{" + std::to_string(ie) + "}";
        latex.DrawLatex(0.22 + 0.07 * ie, 0.92, result.c_str());
      }

      // eta_hist[ieta]->Rebin(5);
      int bin = eta_hist[ieta]->FindBin(0.2);
      if (eta_hist[ieta]->GetBinContent(bin) == 0)
      {
        ieta++;
        std::cout << "no data for eta=" << ieta << std::endl;
        continue;
      }
      eta_hist[ieta]->Scale(1.0 / eta_hist[ieta]->GetBinContent(bin) * pow(10, ie));
      eta_hist[ieta]->Draw("same hist");
      eta_hist[ieta]->SetLineColor(kRed);

      ieta++;
    }

    gPad->SetLogy();

    latex.SetTextSize(0.04);

    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.6, 0.87, "#it{#bf{sPHENIX}} Internal");
    latex.SetTextColor(kRed);
    latex.DrawLatex(0.6, 0.83, "eta integrated");
    latex.SetTextColor(kBlue);
    latex.DrawLatex(0.6, 0.79, "failed QA");

    c1->SaveAs((path + "/tsc_spec" + std::to_string(iplot) + ".pdf").c_str());
    delete c1;
  }

  TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
  h_fail->Draw("COLZ");

  latex.SetTextColor(kBlack);
  std::string stuff = std::to_string(h_fail->GetEntries()) + " failed QA";
  latex.DrawLatex(0.2, 0.87, stuff.c_str());

  c2->SaveAs((path + "/failQA.pdf").c_str());

  TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
  h_expCalib->Draw("COLZ");
  h_expCalib->GetZaxis()->SetRangeUser(0.5, 2);
  gPad->SetRightMargin(0.15);

  latex.SetTextColor(kBlack);

  c4->SaveAs((path + "/expCalib.pdf").c_str());

  h_hot->Write();
  h_hits->Write();
  h_hits_clean->Write();
  h_corrPat->Write();
  h_expCalib->Write();

  f_temp->Close();
}

void LiteCaloEval::draw_spectra(const char *outfile)
{
 
 if (calotype == LiteCaloEval::NONE)
  {
    std::cout << "Did not enter correct calotype. Exiting macro..." << std::endl;
    exit(-1);
  }

  TFile *fout = new TFile(outfile,"UPDATE");

  TH1F *h = nullptr;
  std::string histName;

  TH1F *h_etaSlice = nullptr;
  std::string etaSliceName;

  TH1F *h_cln = nullptr; //for tower clone
  std::string h_cln_nm;

  TH1F *h_es_cln = nullptr; //for es clone
  std::string h_es_nm;


  float scale = 1.0;
  float otherscale = 1.0;

  float fitMin = getFitMin();
  float fitMax = getFitMax();

  double xaxisRange = 2.0;

  if (calotype == LiteCaloEval::HCALIN)
  {
    xaxisRange = 1.0;
  }

  if (calotype == LiteCaloEval::HCALIN || calotype == LiteCaloEval::HCALOUT)
  {
    std::cout << "Drawing hcal tower spectra" << std::endl;

    int starteta = 0;
    int maxeta = 12;
    int power = 23;

    TCanvas *c = new TCanvas();
    c->Divide(2);  // split canvas in two. Left pad holds ieta 0 - 11. Right holds 12-23
    c->SetName("hcal_spectra");

    int cntr = 1;  // allows histos to be separated by 10^2

    // this for loop places you in left or right pad
    for (int k = 1; k < 3; k++)
    {
      c->cd(k);
      gPad->SetLogy(1);

      if (k == 2)
      {
        power = 47;
        starteta = 12;
        maxeta = 24;
      }

      TLegend *t = new TLegend(0.5, 0.8, 0.75, 0.9);
      t->AddEntry("", (boost::format("ieta %d - %d") % starteta % (maxeta - 1)).str().c_str());

      // eta loop
      for (int i = starteta; i < maxeta; i++)
      {
        // get eta slice histogram and overlay in red onto tower spectra
        if (calotype == LiteCaloEval::HCALOUT)
        {
          etaSliceName = "hcalout_eta_" + std::to_string(i);
        }
        else if (calotype == LiteCaloEval::HCALIN)
        {
          etaSliceName = "hcalin_eta_" + std::to_string(i);
        }

        h_etaSlice = (TH1F *) fout->Get(etaSliceName.c_str());

        if (!h_etaSlice)
        {
          std::cout << "ERROR! Could not get hcal eta slice histogram " << i << " ." << std::endl;
          gSystem->Exit(1);
        }

        if (h_etaSlice->GetEntries() == 0.0)
        {
          std::cout << "WARNING! hcal eta slice " << i << " has no entries!" << std::endl;
	  continue;
        }

	h_es_cln = (TH1F *)h_etaSlice->Clone();

	h_es_nm = "cln_es_" + std::to_string(i);

	h_es_cln->SetName(h_es_nm.c_str());


        // phi loop
        for (int j = 0; j < 64; j++)
        {
          if (calotype == LiteCaloEval::HCALOUT)
          {
            histName = "hcal_out_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);
          }
          else if (calotype == LiteCaloEval::HCALIN)
          {
            histName = "hcal_in_eta_" + std::to_string(i) + "_phi_" + std::to_string(j);
          }

          h = (TH1F *) fout->Get(histName.c_str());

          if (!h)
          {
            std::cout << "ERROR! Could not find tower " << histName << "." << std::endl;
            gSystem->Exit(1);
          }

          if (h->GetEntries() == 0.0)
          {
            std::cout << "WARNING! No entries in hcal (" << i << "," << j << ").  Skipping this tower" << std::endl;
            continue;
          }

	  h_cln = (TH1F *)h->Clone();

	  h_cln_nm = (boost::format("h_cln_eta_%d_phi_%d") % i % j).str();

	  h_cln->SetName(h_cln_nm.c_str());


          if (i == 0)
          {
            scale = pow(10, power - i);
          }
          else
          {
            scale = pow(10, power - i - cntr);
          }

          otherscale = h_cln->GetBinContent(1);
	  
          h_cln->Scale(h_es_cln->Integral(h_es_cln->FindBin(fitMin), h_es_cln->FindBin(fitMax) ) / h_cln->Integral(h_cln->FindBin(fitMin), h_cln->FindBin(fitMax) ));
          h_cln->Scale(scale / otherscale);
          h_cln->GetXaxis()->SetRangeUser(0., xaxisRange);
          h_cln->GetYaxis()->SetRangeUser(0.00001, 1e27);
          h_cln->Draw("same hist");

          h_cln = nullptr;

        }  // phi loop

        // draw the combined eta slice spectrum

        h_es_cln->Scale(scale / otherscale);
        h_es_cln->GetXaxis()->SetRangeUser(0., xaxisRange);
        h_es_cln->GetYaxis()->SetRangeUser(0.00001, 1e27);
        h_es_cln->SetLineColor(2);
        h_es_cln->SetLineWidth(3);
        h_es_cln->Draw("same hist");

        h_es_cln = nullptr;

        if (i > 0)
        {
          cntr++;
        }

      }  // eta loop

      t->Draw("same");
      t = nullptr;
      delete t;

    }  // k loop

    c->Write();
    c = nullptr;
    delete c;

  }  // end hcal caloflag

  if (calotype == LiteCaloEval::CEMC)
  {
    std::cout << "Drawing emcal tower spectra" << std::endl;

    int starteta = 0;
    int maxeta = 12;

    for (int k = 1; k < 9; k++)
    {
      int power = 36;
      int cntr = 3;  // cntr == 2 allows histos to be separated by 10^2, cntr == 3 gives separation by 10^3, etc.

      if (k > 1)
      {
        starteta = maxeta;
        maxeta = maxeta + 12;
      }

      TCanvas *c = new TCanvas();
      c->SetName((boost::format("emcal_eta%d_%d") % starteta % (maxeta - 1)).str().c_str());
      gPad->SetLogy(1);

      TLegend *t = new TLegend(0.5, 0.8, 0.75, 0.9);
      t->AddEntry("", (boost::format("ieta %d - %d") % starteta % (maxeta - 1)).str().c_str());

      // eta loop
      for (int i = starteta; i < maxeta; i++)
      {
        etaSliceName = "eta_" + std::to_string(i);

        h_etaSlice = (TH1F *) fout->Get(etaSliceName.c_str());

        if (!h_etaSlice)
        {
          std::cout << "ERROR! Could not get emcal eta slice " << i << ". Exiting analysis." << std::endl;
          exit(-1);
        }

        if (h_etaSlice->GetEntries() == 0.0)
        {
          std::cout << "WARNING! EMCal eta slice " << i << " has no entries!" << std::endl;
	  continue;
        }

	h_es_cln = (TH1F *)h_etaSlice->Clone();

	h_es_nm = "cln_es_" + std::to_string(i);

	h_es_cln->SetName(h_es_nm.c_str());

        // phi loop
        for (int j = 0; j < 256; j++)
        {
          histName = "emc_ieta" + std::to_string(i) + "_phi" + std::to_string(j);

          h = (TH1F *) fout->Get(histName.c_str());

          if (!h)
          {
            std::cout << "ERROR! Could not find " << histName << ". Exiting macro." << std::endl;
            exit(-1);
          }
          if (h->GetEntries() == 0.0)
          {
            std::cout << "WARNING! No entries in emcal (" << i << "," << j << ").  Skipping this tower" << std::endl;
            continue;
          }

	  h_cln = (TH1F *)h->Clone();

	  h_cln_nm = (boost::format("h_cln_eta_%d_phi_%d") % i % j).str();

	  h_cln->SetName(h_cln_nm.c_str());

          if (i == 0 || i % 12 == 0)
          {
            scale = pow(10, power);
          }
          else
          {
            scale = pow(10, power - cntr);
          }

          otherscale = h_cln->GetMaximum();

          h_cln->Scale(h_es_cln->Integral(h_es_cln->FindBin(fitMin), h_es_cln->FindBin(fitMax) ) / h_cln->Integral(h_cln->FindBin(fitMin), h_cln->FindBin(fitMax)));
          h_cln->Scale(scale / otherscale);
          h_cln->GetXaxis()->SetRangeUser(0., xaxisRange);
          h_cln->GetYaxis()->SetRangeUser(0.00001, 1e40);
          h_cln->Draw("same hist");

	  h_cln = nullptr;

        }// phi loop

        h_es_cln->Scale(scale / otherscale);
        h_es_cln->GetXaxis()->SetRangeUser(0., xaxisRange);
        h_es_cln->GetYaxis()->SetRangeUser(0.00001, 1e40);
        h_es_cln->SetLineColor(2);
        h_es_cln->SetLineWidth(3);
        h_es_cln->Draw("same hist");

	h_es_cln = nullptr;

        if (!(i % 12 == 0))
        {
          cntr += 3;  // cntr++ gives separation by 100, +=2 gives 1000
        }

      }// eta loop

      t->Draw("same");

      c->Write();

      t = nullptr;
      delete t;

      c = nullptr;
      delete c;

      std::cout << "Drew canvas " << k << std::endl;
    }  // k loop

  }  // end emcal flag


  fout->Close();
  fout = nullptr;
  delete fout;
  
  std::cout << "Drawing histos is complete." << std::endl;

}  // end draw spectra f'n

void LiteCaloEval::fit_info(const char * outfile, const int runNum)
{
  TFile *fout = new TFile(outfile,"UPDATE");

  int eta;
  int phi;
  double towers;
//  int badTowers = 0;

  if (calotype == LiteCaloEval::HCALIN || calotype == LiteCaloEval::HCALOUT)
  {
    eta = 24;
    phi = 64;
    towers = 1536.0;
  }
  else if (calotype == LiteCaloEval::CEMC)
  {
    eta = 96;
    phi = 256;
    towers = 24576.0;
  }
  else
  {
    std::cout << "calotype not set. Exiting." << std::endl;
    exit(-1);
  }

  TH2F *errMap = new TH2F("errMap", "", eta, 0, eta, phi, 0, phi);
  errMap->GetXaxis()->SetTitle("#eta Bin");
  errMap->GetYaxis()->SetTitle("#phi Bin");

  // make chi squared/ndf plots
  TH1F *chi2 = new TH1F("chi2", "", 500, 0, 50);
  chi2->GetXaxis()->SetTitle("#chi^{2} / NDF");
  chi2->GetYaxis()->SetTitle("Counts");

  // make chi squared/ndf map
  TH2F *chi2Map = new TH2F("chi2Map", "", eta, 0, eta, phi, 0, phi);
  chi2Map->GetXaxis()->SetTitle("#eta Bin");
  chi2Map->GetYaxis()->SetTitle("#phi Bin");

  // map of tower with failed fits
  TH2F *fitFail = new TH2F("fitFail", "", eta, 0, eta, phi, 0, phi);
  fitFail->GetXaxis()->SetTitle("#eta bin");
  fitFail->GetYaxis()->SetTitle("#phi bin");

  TString histname;
  TH1F *htmp = nullptr;
  TF1 *fn = nullptr;

  //testing stuff
  TH2F *cp = nullptr; // to hold corrpat 
  cp = (TH2F *)fout->Get("corrPat");

  if(!cp)
    {
      std::cout << "Error! Did not get corrPat histogram. Exiting analysis." << std::endl;
      exit(-1);
    }

  TGraphErrors *tmp_avgTSC = new TGraphErrors();
  tmp_avgTSC->SetName("g_avgTSC");
  tmp_avgTSC->GetXaxis()->SetTitle("run number");
  tmp_avgTSC->GetYaxis()->SetTitle("Mean Towerslope Correction"); 
  tmp_avgTSC->SetMarkerStyle(8);
  tmp_avgTSC->SetMarkerSize(1);
    

  double sum4avg = 0.0;
  double tscAvg = 0.0;
  double sum4SE = 0.0;
  double SE = 0.0;

  // phi loop
  for (int i = 0; i < eta; i++)
    {
      // eta loop
      for (int j = 0; j < phi; j++)
	{
	  // for ohcal
	  if (calotype == LiteCaloEval::HCALOUT)
	    {
	      histname = (boost::format("hcal_out_eta_%d_phi_%d") % i % j).str();
	    }

	  // ihcal
	  if (calotype == LiteCaloEval::HCALIN)
	    {
	      histname = (boost::format("hcal_in_eta_%d_phi_%d") % i % j).str();
	    }

	  // emcal
	  if (calotype == LiteCaloEval::CEMC)
	    {
	      histname = (boost::format("emc_ieta%d_phi%d") % i % j).str();
	    }

	  sum4avg += (cp->GetBinContent(i+1, j+1));

	  htmp = (TH1F *) fout->Get(histname.Data());

	  fn = (TF1 *) htmp->GetFunction("myexpo");

	  errMap->SetBinContent(i + 1, j + 1, fn->GetParError(1));

	  chi2->Fill(fn->GetChisquare() / fn->GetNDF());

	  chi2Map->SetBinContent(i + 1, j + 1, fn->GetChisquare() / fn->GetNDF());

	  if (fn->GetChisquare() / fn->GetNDF() > 5)
	    {
	      fitFail->Fill(i, j);
//	      badTowers++;
	    }

	}  // end inner forloop

    }  // end outer forloop

  tscAvg = sum4avg / towers;

  //get standard error of avg TSC
  for(int i = 0; i < eta; i++)
    {
      for(int j = 0; j < phi; j++)
	{
	  sum4SE += pow(cp->GetBinContent(i+1, j+1) - tscAvg, 2); //getting part of the std dev
	}
    }

  SE = sqrt(sum4SE/towers) / sqrt(towers);

  tmp_avgTSC->SetPoint(0, runNum, tscAvg); 
  tmp_avgTSC->SetPointError(0, 0.0, SE);

  errMap->Write();
  chi2->Write();
  chi2Map->Write();
  fitFail->Write();
  tmp_avgTSC->Write();

  errMap = nullptr;
  delete errMap;

  chi2 = nullptr;
  delete chi2;

  chi2Map = nullptr;
  delete chi2Map;

  fitFail = nullptr;
  delete fitFail;

  fout->Close();
  fout = nullptr;
  delete fout;

  std::cout << "Finished fit info" << std::endl;
}
