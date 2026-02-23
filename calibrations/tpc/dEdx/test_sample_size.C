#include "GlobaldEdxFitter.h"

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TFile.h>

#include <algorithm>
#include <numeric>

void test_sample_size(const std::string& infile="/sphenix/tg/tg01/hf/mjpeters/run53877_tracks/track_output_53877_*.root")
{
  std::vector<float> samplesizes = {1000,2000,5000,10000,20000};//,50000,100000,200000,500000,1000000};//,2000000,5000000};

  const int n_samples = 20;
  const float fluctuation_ymin = 5.;
  const float fluctuation_ymax = 26.;
  const int distribution_nbins = 30;
  const float distribution_xmin = 5.;
  const float distribution_xmax = 26.;

  std::vector<std::vector<float>> fitvalues_all;
  std::vector<float> fitvalues_avg;
  std::vector<float> fitvalues_stdev;

  std::vector<TCanvas*> fluctuations;
  std::vector<TCanvas*> distributions;

  std::vector<TH1F*> dist_h;

  std::vector<std::unique_ptr<GlobaldEdxFitter>> gfs;

  for(int i=0; i<samplesizes.size(); i++)
  {
    fitvalues_all.emplace_back();
    gfs.push_back(std::make_unique<GlobaldEdxFitter>());

    const std::string& fluctuation_canvasname = "fluctuations_"+std::to_string((int)floor(samplesizes[i]));
    const std::string& distribution_canvasname = "distributions_"+std::to_string((int)floor(samplesizes[i]));
    fluctuations.push_back(new TCanvas(fluctuation_canvasname.c_str(),fluctuation_canvasname.c_str(),600,600));
    distributions.push_back(new TCanvas(distribution_canvasname.c_str(),distribution_canvasname.c_str(),600,600));

    for(int j=0;j<n_samples;j++)
    {
      gfs[i]->processResidualData(infile,floor(samplesizes[i]),j*samplesizes[i]);
      double min = gfs[i]->get_minimum();
      std::cout << "minimum: " << min << std::endl;
      fitvalues_all[i].push_back(min);
      if(j<n_samples-1)
      {
        gfs[i]->reset();
      }
/*
      tf1s[i]->cd();
      TF1* tf1copy = gfs[i]->create_TF1(("ntrk_"+std::to_string(samplesizes[i])).c_str());
      tf1copy->SetLineColor(base_color);
      tf1copy->GetYaxis()->SetRangeUser(1.,tf1copy->GetMaximum());
      if(i==0) tf1copy->Draw();
      else tf1copy->Draw("SAME");
*/
    }
  }

  std::vector<float> sample_index(n_samples);
  std::iota(sample_index.begin(),sample_index.end(),0.);

  for(int i=0; i<samplesizes.size(); i++)
  {
    fluctuations[i]->cd();
    TGraph* g = new TGraph(n_samples,sample_index.data(),fitvalues_all[i].data());
    g->GetYaxis()->SetRangeUser(fluctuation_ymin,fluctuation_ymax);
    g->SetMarkerStyle(kFullCircle);
    g->SetMarkerSize(1.);
    g->Draw("APL");

    distributions[i]->cd();
    std::string hname = "h_"+std::to_string(floor(samplesizes[i]));
    std::string htitle = "Distribution of fit results for sample size "+std::to_string(floor(samplesizes[i]));
/*
    auto bounds = std::minmax_element(fitvalues_all[i].begin(),fitvalues_all[i].end());
    float lowerbound = floor(*bounds.first);
    float upperbound = ceil(*bounds.second);
*/
    TH1F* h = new TH1F(hname.c_str(),htitle.c_str(),distribution_nbins,distribution_xmin,distribution_xmax);
    for(int j=0; j<n_samples; j++)
    {
      h->Fill(fitvalues_all[i][j]);
    }
    h->Draw();
  }

  for(int i=0; i<samplesizes.size(); i++)
  {
    float avg = 0.;
    float stdev = 0.;
    for(int j=0; j<n_samples; j++)
    {
      avg += fitvalues_all[i][j];
    }
    avg /= (float)n_samples;

    for(int j=0; j<n_samples; j++)
    {
      stdev += pow(fitvalues_all[i][j]-avg,2.);
    }
    stdev /= (float)n_samples;
    stdev = std::sqrt(stdev);

    fitvalues_avg.push_back(avg);
    fitvalues_stdev.push_back(stdev);
  }

  std::vector<float> errx(n_samples,0.);

  TCanvas* cg = new TCanvas("cg","sizes",600,600);
  TGraph* g = new TGraphErrors(samplesizes.size(),samplesizes.data(),fitvalues_avg.data(),errx.data(),fitvalues_stdev.data());
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerSize(1);
  g->Draw("APL");
  cg->SetLogx();

  TCanvas* cbg = new TCanvas("vsbetagamma","vsbetagamma",600,600);
  TGraph* gbg = gfs.back()->graph_vsbetagamma(fitvalues_avg.back());
  gbg->SetMarkerStyle(kFullCircle);
  gbg->SetMarkerSize(0.2);
  gbg->Draw("AP");
  cbg->SetLogx();

  double best_A = fitvalues_avg.back();

  TF1* bethe = new TF1("bethebloch_vslnbg",bethe_bloch_new_1D_wrapper,0.,100.,2,1);
  bethe->SetParameter(0,best_A);
  bethe->SetNpx(1000);
  bethe->Draw("SAME");

  TF1* bethe_directfit = new TF1("bethebloch_directfit",bethe_bloch_new_1D_wrapper,0.,10.,1,1);
  bethe_directfit->SetParameter(0,best_A);
  bethe_directfit->SetLineColor(kBlue);
  gbg->Fit(bethe_directfit);
  double newbest_A = bethe_directfit->GetParameter(0);
  std::cout << "new best: " << newbest_A << std::endl;

  TCanvas* cbands = new TCanvas("bands","bands",600,600);
  TGraph* gp = gfs.back()->graph_vsp();
  gp->SetMarkerStyle(kFullCircle);
  gp->SetMarkerSize(0.1);
  gp->Draw("AP");
  cbands->SetLogx();

  for(double mass : {dedx_constants::m_pi, dedx_constants::m_K, dedx_constants::m_p, dedx_constants::m_d})
  {
    TF1* band = new TF1(("band_"+std::to_string(mass)).c_str(),bethe_bloch_vs_p_wrapper_new_1D,0.,10.,2,1);
    band->SetParameter(0,best_A);
    band->SetParameter(1,mass);
    band->SetNpx(1000);
    band->Draw("SAME");

    TF1* directband = new TF1(("directband_"+std::to_string(mass)).c_str(),bethe_bloch_vs_p_wrapper_new_1D,0.,10.,2,1);
    directband->SetLineColor(kBlue);
    directband->SetParameters(best_A,mass);
    directband->SetNpx(1000);
    directband->Draw("SAME");
  }

  TCanvas* cb = new TCanvas("fullbands","fullbands",600,600);
  TFile* f_h = TFile::Open("/sphenix/tg/tg01/hf/mjpeters/run53877_tracks/dedx/merged_dedx.root");
  TH2F* dedx_h = (TH2F*)f_h->Get("dedx_log_30");
  dedx_h->Draw("COLZ");
  cb->SetLogz();

  for(double mass : {dedx_constants::m_pi, dedx_constants::m_K, dedx_constants::m_p, dedx_constants::m_d})
  {
    TF1* band = new TF1(("band_"+std::to_string(mass)).c_str(),bethe_bloch_vs_logp_wrapper_new_1D,-1.,5.,2,1);
    band->SetParameter(0,best_A);
    band->SetParameter(1,mass);
    band->Draw("SAME");

    TF1* directband = new TF1(("directband_"+std::to_string(mass)).c_str(),bethe_bloch_vs_logp_wrapper_new_1D,-1.,5.,2,1);
    directband->SetLineColor(kBlue);
    directband->SetParameters(newbest_A,mass);
    directband->SetNpx(1000);
    directband->Draw("SAME");
  }

  TFile* fout = new TFile("dedxfitvals.root","RECREATE");
  for(auto& c : fluctuations)
  {
    c->Write();
  }
  for(auto& c : distributions)
  {
    c->Write();
  }
  cg->Write();
  cbg->Write();
  cbands->Write();
  cb->Write();
  fout->Close();
}
