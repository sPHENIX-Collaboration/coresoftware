#include <cdbobjects/CDBTTree.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TStyle.h>
#include <TText.h>
#include <TTree.h>

#include <boost/format.hpp>

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//////////////
// Functions //
//////////////
double SingleGaussianFit(TH1D *hist, double &mean1, double &sigma1);
// Function to do fitting return value : value for HotChannel cut

/////////////////////
// Global parameters//
/////////////////////
// int chip = 26;
int mod = 14;
double sig_cut = 3.0;
bool Writecsv = false;
bool debug = false;
// chip : # of chips on each half ladder(defalt : 26)
// mod : # of moudles for each felix server(default : 14)
// sig_cut : sigma cut for hot/cold channel determination default : 3.0
// Writecsv : flag to write csv file (default : false )
//  debug : used to turn ON/OFF cout
/////////////
// File Path//
/////////////
// std::string map_input_path =  "/sphenix/tg/tg01/commissioning/INTT/QA/hitmap/2024/";
// std::string cdb_output_path = "/sphenix/tg/tg01/commissioning/INTT/QA/hotdeadmap/CDB/2024/";
// std::string root_output_path ="/sphenix/tg/tg01/commissioning/INTT/QA/hotdeadmap/rootfile/2024/";
// std::string csv_output_path = "/sphenix/tg/tg01/commissioning/INTT/QA/hotdeadmap/csv_file/2024/";
std::string map_input_path = "./";
std::string cdb_output_path = "./";
std::string root_output_path = "./";
std::string csv_output_path = "./";
// map_input_path : location of the HitMap file created by InttHitMap.cc/h
// cdb_output_path : cdb output file path
// root_output_path : output file path for QA
// csv_output_path : csv output file path (used for Grafana online monitoring)

struct Half_Chip
{
  int felix_id_;
  int module_id_;
  int chip_id_;
};

std::vector<Half_Chip> half_chips =  // Chip number 0~25
    {
        // Felix 0
        {0, 7, 14},
        // Felix 2
        {2, 9, 15},
        // Felix 3
        {3, 13, 20},
        {3, 13, 22},
        {3, 13, 24},
        // Felix 7
        {7, 0, 0},
        {7, 0, 1},
        {7, 0, 2},
        {7, 0, 3},
        {7, 0, 4},
        {7, 0, 13},
        {7, 0, 14},
        {7, 0, 15},
        {7, 0, 16},
        {7, 0, 17},
        {7, 0, 6},
        {7, 0, 8},
        {7, 0, 10},
        {7, 0, 12},
        {7, 0, 19},
        {7, 0, 21},
        {7, 0, 23},
        {7, 0, 25},
        {7, 1, 0},
        {7, 1, 1},
};

////////////////////////////
// List of Half entry chips//
////////////////////////////

void InttChannelClassifier(int runnumber = 20869)  // runnumber
{
  gStyle->SetOptFit();
  ////////////////////////////////////////
  // Load HitMap                        //
  ////////////////////////////////////////
  int in_sig = 10 * sig_cut;
  std::string rootfilename = map_input_path + "hitmap_run" + std::to_string(runnumber) + ".root";
  std::string cdbttree_name = cdb_output_path + "cdb_inttbadchanmap_" + std::to_string(runnumber) + ".root";
  TFile *file = nullptr;
  file = TFile::Open(rootfilename.c_str(), "READ");

  ///////////////////////////////////////////////////
  // Create Histogram for Hot and Dead channel map //
  ///////////////////////////////////////////////////

  TH2D *h2_AllMap[8][14];   // Original HitMap
  TH2D *h2_ColdMap[8][14];  // 2D histogram for coldmap
  TH2D *h2_HalfMap[8][14];  // 2D histogram for halfmap (half entry chips)
  TH2D *h2_DeadMap[8][14];  // 2D histogram for deadmap
  TH2D *h2_HotMap[8][14];   // 2D histogram for hotmap
  for (int i = 0; i < 8; i++)
  {
    for (int j = 0; j < 14; j++)
    {
      h2_AllMap[i][j] = (TH2D *) file->Get((boost::format("HitMap_%d_%d") % i % j).str().c_str());
      h2_ColdMap[i][j] = new TH2D((boost::format("ColdMap_%d_%d") % i % j).str().c_str(), (boost::format("ColdMap_%d_%d") % i % j).str().c_str(), 128, 0, 128, 26, 0, 26);
      h2_HalfMap[i][j] = new TH2D((boost::format("HalfMap_%d_%d") % i % j).str().c_str(), (boost::format("HalfMap_%d_%d") % i % j).str().c_str(), 128, 0, 128, 26, 0, 26);
      h2_HotMap[i][j] = new TH2D((boost::format("HotMap_%d_%d") % i % j).str().c_str(), (boost::format("HotMap_%d_%d") % i % j).str().c_str(), 128, 0, 128, 26, 0, 26);
      h2_DeadMap[i][j] = new TH2D((boost::format("DeadMap_%d_%d") % i % j).str().c_str(), (boost::format("DeadMap_%d_%d") % i % j).str().c_str(), 128, 0, 128, 26, 0, 26);
    }
  }
  //////////////////////////////////////////
  // Define condition for the hot channel //
  //////////////////////////////////////////
  double HotChannelCut_A_Fit[8][14] = {{0.}};
  double HotChannelCut_B_Fit[8][14] = {{0.}};
  double ColdChannelCut_A_Fit[8][14] = {{0.}};
  double ColdChannelCut_B_Fit[8][14] = {{0.}};
  double par_meanA[8][14] = {{0.}};   // Array to save the mean & sigma value, [0][module] = mean, [1][module] = sigma
  double par_sigmaA[8][14] = {{0.}};  // Array to save the mean & sigma value, [0][module] = mean, [1][module] = sigma
  double par_meanB[8][14] = {{0.}};   // Array to save the mean & sigma value, [0][module] = mean, [1][module] = sigma
  double par_sigmaB[8][14] = {{0.}};  // Array to save the mean & sigma value, [0][module] = mean, [1][module] = sigma
  /////////////////////////////////////////////////
  // Create TFile and TTree to save information  //
  // These are used to check fitting result      //
  /////////////////////////////////////////////////
  std::string csvfilename = csv_output_path + "NumOfHot.csv";
  std::ofstream csvFile(csvfilename, std::ios::app);
  int NumOfHot = 0;
  if (!csvFile.is_open())
  {
    std::cout << csvfilename << std::endl;
    std::cerr << "Unable to open the file." << std::endl;
    return;
  }
  std::string outputfile = root_output_path + "InttHotDeadMap_test_" + std::to_string(runnumber) + "_" + std::to_string(in_sig) + ".root";
  TFile *sfile = new TFile(outputfile.c_str(), "RECREATE");
  TTree *st = new TTree("tree", "tree");
  double ch_entry, mean_gaus, sigma_gaus = 0.;
  int felix_, module_, chip_id_, chan_, type_, flag_;
  st->Branch("felix", &felix_);
  st->Branch("module", &module_);
  st->Branch("chip_id", &chip_id_);
  st->Branch("chan", &chan_);
  st->Branch("flag", &flag_);
  st->Branch("ch_entry", &ch_entry);
  st->Branch("type", &type_);
  st->Branch("mean", &mean_gaus);
  st->Branch("sigma", &sigma_gaus);

  ///////////////////////////////////
  // Fill the Dead & Cold & HotMap //
  ///////////////////////////////////
  TCanvas *canA[8];
  TCanvas *canB[8];
  for (int i = 0; i < 8; i++)
  {
    canA[i] = new TCanvas((boost::format("TypeA_Felix_%d") % i).str().c_str(), (boost::format("TypeA_Felix_%d") % i).str().c_str(), 1200, 1200);
    canB[i] = new TCanvas((boost::format("TypeB_Felix_%d") % i).str().c_str(), (boost::format("TypeB_Felix_%d") % i).str().c_str(), 1200, 1200);
    canA[i]->Divide(7, 2);
    canB[i]->Divide(7, 2);
  }
  TH1D *h1_hist_fit_A[8][14];
  TH1D *h1_hist_fit_B[8][14];
  double mean_first = 0;
  for (int felix = 0; felix < 8; felix++)
  {
    for (int i = 0; i < 14; i++)
    {
      h1_hist_fit_A[felix][i] = new TH1D((boost::format("h1_hist_fit_A%d_%d") % felix % i).str().c_str(), (boost::format("h1_hist_fit_A%d_%d") % felix % i).str().c_str(), 100, 0, 0.03);
      h1_hist_fit_B[felix][i] = new TH1D((boost::format("h1_hist_fit_B%d_%d") % felix % i).str().c_str(), (boost::format("h1_hist_fit_B%d_%d") % felix % i).str().c_str(), 100, 0, 0.03);
      for (int j = 0; j < 26; j++)
      {
        for (int chan = 0; chan < 128; chan++)
        {
          if (j < 5 || (j > 12 && j < 18))  // Condition for type B
          {
            h1_hist_fit_B[felix][i]->Fill(h2_AllMap[felix][i]->GetBinContent(chan + 1, j + 1));
          }
          else
          {
            h1_hist_fit_A[felix][i]->Fill(h2_AllMap[felix][i]->GetBinContent(chan + 1, j + 1));
          }
        }
      }
      double mean, sigma;
      canA[felix]->cd(i + 1);
      HotChannelCut_A_Fit[felix][i] = SingleGaussianFit(h1_hist_fit_A[felix][i], mean, sigma);
      ColdChannelCut_A_Fit[felix][i] = mean - sig_cut * sigma;
      if (mean_first == 0)
      {
        mean_first = mean;
      }
      if (mean < 0.005)
      {
        HotChannelCut_A_Fit[felix][i] = 1;
        ColdChannelCut_A_Fit[felix][i] = 1;
        if (mean - sig_cut * sigma < 0)
        {
          mean = -1;
          sigma = -1;
        }
      }
      par_meanA[felix][i] = mean;
      par_sigmaA[felix][i] = sigma;
      canB[felix]->cd(i + 1);
      HotChannelCut_B_Fit[felix][i] = SingleGaussianFit(h1_hist_fit_B[felix][i], mean, sigma);
      ColdChannelCut_B_Fit[felix][i] = mean - sig_cut * sigma;
      if (mean < 0.005)
      {
        HotChannelCut_B_Fit[felix][i] = 1;
        ColdChannelCut_B_Fit[felix][i] = 1;
        if (mean - sig_cut * sigma < 0)
        {
          mean = -1;
          sigma = -1;
        }
      }
      par_meanB[felix][i] = mean;
      par_sigmaB[felix][i] = sigma;
    }
  }

  CDBTTree *cdbttree = new CDBTTree(cdbttree_name);
  int size = 0;
  sfile->cd();
  TDirectory *dir[8];
  for (int felix = 0; felix < 8; felix++)
  {
    dir[felix] = sfile->mkdir((boost::format("Felix_%d") % felix).str().c_str());
    dir[felix]->cd();
    for (int i = 0; i < 14; i++)
    {
      if (debug)
      {
        std::cout << "moudle : " << i << " Type A " << HotChannelCut_A_Fit[felix][i] << std::endl;
        std::cout << "moudle : " << i << " Type B " << HotChannelCut_B_Fit[felix][i] << std::endl;
      }
      for (int j = 0; j < 26; j++)
      {
        if (debug)
        {
          std::cout << "Felix : " << felix << " moudle : " << i << " Type A and chip : " << j << "  " << HotChannelCut_A_Fit[felix][i] << std::endl;
        }
        for (int chan = 0; chan < 128; chan++)
        {
          // double entry = h1_chip[i][j]->GetBinContent(chan + 1);
          double entry = h2_AllMap[felix][i]->GetBinContent(chan + 1, j + 1);
          felix_ = felix;
          ch_entry = entry;
          module_ = i;
          chip_id_ = j;
          chan_ = chan;
          flag_ = 0;
          // flag = 0 : Good channels
          // flag = 1 : Dead channels
          // flag = 2 : Half channels
          // flag = 4 : Cold channels
          // flag = 8 : Hot channels
          if (j < 5 || (j > 12 && j < 18))  // Condition for type B
          {
            type_ = 1;
            mean_gaus = par_meanB[felix][i];
            sigma_gaus = par_sigmaB[felix][i];
          }
          else
          {
            type_ = 0;
            mean_gaus = par_meanA[felix][i];
            sigma_gaus = par_sigmaA[felix][i];
          }
          ///////////////////////////////////////////////////////////
          //                 Half entry chips                      //
          ///////////////////////////////////////////////////////////
          bool is_half = false;
          for (const auto &chip : half_chips)
          {
            if (chip.felix_id_ == felix && chip.module_id_ == i && chip.chip_id_ == j)
            {
              is_half = true;
              break;
            }
          }
          if (is_half)
          {
            h2_HalfMap[felix][i]->Fill(chan, j);
            flag_ = 2;
          }
          ///////////////////////////////////////////////////////////
          //                 Dead Channel selection                //
          ///////////////////////////////////////////////////////////
          else if (entry == 0)
          {
            h2_DeadMap[felix][i]->Fill(chan, j);
            flag_ = 1;
          }
          else
          {
            ///////////////////////////////////////////////////////////
            //           Hot Cold Channel selection For tpye B       //
            ///////////////////////////////////////////////////////////
            if (j < 5 || (j > 12 && j < 18))  // Condition for type B
            {
              if (entry > HotChannelCut_B_Fit[felix][i])
              {
                h2_HotMap[felix][i]->Fill(chan, j);
                flag_ = 8;
                NumOfHot++;
              }
              if (entry < ColdChannelCut_B_Fit[felix][i])
              {
                h2_ColdMap[felix][i]->Fill(chan, j);
                flag_ = 4;
              }
            }
            else  // For tpye A
            {
              ///////////////////////////////////////////////////////////
              //            Hot Cold Channel selection for tpye A      //
              ///////////////////////////////////////////////////////////
              if (entry > HotChannelCut_A_Fit[felix][i])
              {
                h2_HotMap[felix][i]->Fill(chan, j);
                flag_ = 8;
                NumOfHot++;
              }
              if (entry < ColdChannelCut_A_Fit[felix][i])
              {
                h2_ColdMap[felix][i]->Fill(chan, j);
                flag_ = 4;
              }
            }
          }
          ////////////////////////////////////////////////////////////////
          //           Fill CDBTTree if channel is not good             //
          ////////////////////////////////////////////////////////////////
          if (flag_ != 0)
          {
            cdbttree->SetIntValue(size, "felix_server", felix);
            cdbttree->SetIntValue(size, "felix_channel", module_);
            cdbttree->SetIntValue(size, "chip", chip_id_);
            cdbttree->SetIntValue(size, "channel", chan_);
            cdbttree->SetIntValue(size, "flag", flag_);
            ++size;
          }
          st->Fill();
        }
      }
      h2_AllMap[felix][i]->Write();
      h2_DeadMap[felix][i]->Write();
      h2_ColdMap[felix][i]->Write();
      h2_HalfMap[felix][i]->Write();
      h2_HotMap[felix][i]->Write();
      h1_hist_fit_A[felix][i]->Write();
      h1_hist_fit_B[felix][i]->Write();
    }
    sfile->cd();
    canA[felix]->Write();
    canB[felix]->Write();
  }

  cdbttree->SetSingleIntValue("size", size);
  cdbttree->Commit();
  cdbttree->CommitSingle();
  cdbttree->WriteCDBTTree();
  st->Write();
  // Add content to the end of the file
  if (Writecsv)
  {
    csvFile << runnumber << "," << NumOfHot << "\n";
  }

  // Close the file
  csvFile.close();

  sfile->Close();
  file->Close();
}

double SingleGaussianFit(TH1D *hist, double &mean1, double &sigma1)
{
  // backup
  // TF1* g1 = new TF1("g1", "gaus",0,0.005);
  // TF1* g2 = new TF1("g2", "gaus",0.004,0.015);
  // TF1* g3 = new TF1("g3", "gaus",0.09,0.02);
  //
  double hist_mean = hist->GetMean();
  double hist_std = hist->GetStdDev();
  TF1 *SingleGaussian[7];
  SingleGaussian[0] = new TF1("singleGaussian0", "gaus", hist_mean - 1.3 * hist_std, hist_mean + 1.3 * hist_std);
  SingleGaussian[1] = new TF1("singleGaussian1", "gaus", hist_mean - 2 * 1.3 * hist_std, 1.3 * hist_mean);
  SingleGaussian[2] = new TF1("singleGaussian2", "gaus", hist_mean, hist_mean + 2 * 1.3 * hist_std);
  // SingleGaussian[3] = new TF1("singleGaussian3", "gaus", 0.012, 0.0017);
  // SingleGaussian[4] = new TF1("singleGaussian4", "gaus", 0.010, 0.015);
  // SingleGaussian[5] = new TF1("singleGaussian5", "gaus", 0.006, 0.008);
  // SingleGaussian[6] = new TF1("singleGaussian6", "gaus", 0.04, 0.055);
  TF1 *FirstGaussian = new TF1("First_Gaussian", "gaus", 0.003, 0.03);
  //  double constant = 0;
  double chi2 = 0;
  //  double sigma0 = 0;
  //  double mean0 = 0;
  int ndf = 0;
  double chi2ndf = 0;
  //  double par[9];
  //  double sigma_max = 1e-3;
  //  double sigma_min = 1e-4;
  bool DoMultifit = true;  // By default, try to do fitting with several Gaussian
  hist->Fit(FirstGaussian, "RQ");
  mean1 = FirstGaussian->GetParameter(1);
  chi2 = FirstGaussian->GetChisquare();
  ndf = FirstGaussian->GetNDF();
  sigma1 = FirstGaussian->GetParameter(2);
  chi2ndf = chi2 / ndf;
  int labbel = -1;
  //  int flag = -1;
  // double _mean[7] = {0};
  // double _constant[7] = {0};
  // double _chi2[7] = {0};
  // double _chi2ndf[7] = {0};
  // double _ndf[7] = {0};
  // double _sigma[7] = {0};
  //  int _flag[7] = {0};
  std::array<double, 7> _mean{};
  std::array<double, 7> _constant{};
  std::array<double, 7> _chi2{};
  std::array<double, 7> _chi2ndf{};
  std::array<double, 7> _sigma{};
  std::array<double, 7> _ndf{};
  std::array<int, 7> _flag{};
  _mean.fill(-1);
  _constant.fill(-1);
  _chi2.fill(-1);
  _chi2ndf.fill(-1);
  _sigma.fill(-1);
  _ndf.fill(-1);
  _flag.fill(-1);

  if (DoMultifit)
  {
    hist->Fit(SingleGaussian[0], "RQ");
    hist->Fit(SingleGaussian[1], "RQ+");
    hist->Fit(SingleGaussian[2], "RQ+");

    mean1 = 0;
    sigma1 = 0;

    for (int i = 0; i < 3; i++)
    {
      _mean.at(i) = SingleGaussian[i]->GetParameter(1);
      _constant.at(i) = SingleGaussian[i]->GetParameter(1);
      _chi2.at(i) = SingleGaussian[i]->GetChisquare();
      _ndf.at(i) = SingleGaussian[i]->GetNDF();
      _chi2ndf.at(i) = _chi2.at(i) / _ndf.at(i);
      _sigma.at(i) = SingleGaussian[i]->GetParameter(2);
      if (_chi2ndf.at(i) > 100)
      {
        _mean.at(i) = -1;
        _flag.at(i) = 0;
      }
      if (_sigma.at(i) < 0.00001)
      {
        _mean.at(i) = -1;
        _flag.at(i) = 1;
      }
      if (_mean.at(i) > 0.025 || SingleGaussian[i]->Eval(_mean.at(i)) < 55 || SingleGaussian[i]->Eval(_mean.at(i)) > 2000)
      {
        _mean.at(i) = -1;
        _flag.at(i) = 2;
      }
    }
    mean1 = _mean.at(0);
    sigma1 = _sigma.at(0);
    //    flag = _flag.at(0);
    for (int i = 1; i < 3; i++)
    {
      if (mean1 < _mean.at(i) && _sigma.at(i) != 0)
      {
        mean1 = _mean.at(i);
        sigma1 = _sigma.at(i);
        chi2ndf = _chi2ndf.at(i);
        labbel = i;
      }
    }
  }

  //  flag = _flag.at(labbel);
  TText *text = new TText(0.7, 0.7, (boost::format("chi2/ndf: %.1f , %d") % chi2ndf % labbel).str().c_str());
  text->SetNDC();
  text->SetTextSize(0.03);
  text->Draw("SAME");

  TText *text2 = new TText(0.7, 0.65, (boost::format("sigma: %.10f") % sigma1).str().c_str());
  text2->SetNDC();
  text2->SetTextSize(0.03);
  text2->Draw("SAME");

  TText *text3 = new TText(0.7, 0.6, (boost::format("mean: %.4f") % mean1).str().c_str());
  text3->SetNDC();
  text3->SetTextSize(0.03);
  text3->Draw("SAME");

  TText *text4 = new TText(0.7, 0.55, (boost::format("sigma: %.6f") % sigma1).str().c_str());
  text4->SetNDC();
  text4->SetTextSize(0.03);
  text4->Draw("SAME");
  /* used for debugging purpose.
    TText *text5 = new TText(0.7, 0.50, (boost::format("%d %d %d %d %d %d %d") % _flag.at[0] %_flag[1] %_flag[2] %_flag[3] %_flag[4] %_flag[5] %_flag[6]).str().c_str());
    text5->SetNDC();
    text5->SetTextSize(0.03);
    text5->Draw("SAME");
    text5->Draw("SAME");
  */
  /////////////////////////////
  // Definiton of Hot Channel //
  /////////////////////////////
  TLine *hotLine = new TLine(mean1 + sigma1 * sig_cut, 0, mean1 + sigma1 * sig_cut, hist->GetMaximum());
  hotLine->SetLineColor(kRed);
  hotLine->SetLineWidth(2);
  hotLine->Draw("SAME");
  /////////////////////////////
  // Definiton of Cold Channel//
  /////////////////////////////
  TLine *coldLine = new TLine(mean1 - sigma1 * sig_cut, 0, mean1 - sigma1 * sig_cut, hist->GetMaximum());
  coldLine->SetLineColor(kBlue);
  coldLine->SetLineWidth(2);
  coldLine->Draw("SAME");
  double cutvalue = mean1 + fabs(sigma1) * sig_cut;
  return cutvalue;  // return value : hot channel cut for this ladders
}
