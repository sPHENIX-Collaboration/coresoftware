#include "HCalCosmics.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include "Math/SpecFuncMathCore.h"

#include <Event/Event.h>
#include <Event/packet.h>
#include <TCanvas.h>
#include <cassert>
#include <sstream>

using namespace std;

HCalCosmics::HCalCosmics(const std::string &name, const std::string &fname)
  : SubsysReco(name)
  , detector("HCALIN")
  , outfilename(fname)
{
}

HCalCosmics::~HCalCosmics()
{
}

int HCalCosmics::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << std::endl
            << "HCalCosmics::Init" << std::endl;
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  for (int ieta = 0; ieta < n_etabin; ++ieta)
  {
    for (int iphi = 0; iphi < n_phibin; ++iphi)
    {
      std::string channel_histname = "h_channel_" + std::to_string(ieta) + "_" + std::to_string(iphi);
      h_channel_hist[ieta][iphi] = new TH1F(channel_histname.c_str(), "", 100, 0, 10000);
    }
  }
  h_waveformchi2 = new TH2F("h_waveformchi2", "", 1000, 0, 10000, 1000, 0, 100000);
  h_waveformchi2->GetXaxis()->SetTitle("peak (ADC)");
  h_waveformchi2->GetYaxis()->SetTitle("chi2");
  h_mip = new TH1F("h_mip", "", 100, 0, 10000);

  h_time_energy = new TH2F("h_time_energy", "", 100, -10, 10, 100, -50, 1e3);


  event = 0;
  return 0;
}

int HCalCosmics::process_event(PHCompositeNode *topNode)
{
  if (event % 10000 == 0)
  {
    std::cout << "HCalCosmics::process_event " << event << std::endl;
  }
  process_towers(topNode);
  event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int HCalCosmics::process_towers(PHCompositeNode *topNode)
{
  ostringstream nodenamev2;
  nodenamev2.str("");
  nodenamev2 << prefix << detector;

  TowerInfoContainer *towers = findNode::getClass<TowerInfoContainer>(topNode, nodenamev2.str());
  if (!towers)
  {
    std::cout << std::endl
              << "Didn't find node " << nodenamev2.str() << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int size = towers->size();
  for (int channel = 0; channel < size; channel++)
  {
    TowerInfo *tower = towers->get_tower_at_channel(channel);
    float energy = tower->get_energy();
    float chi2 = tower->get_chi2();
    float time = tower->get_time_float();
    unsigned int towerkey = towers->encode_key(channel);
    int ieta = towers->getTowerEtaBin(towerkey);
    int iphi = towers->getTowerPhiBin(towerkey);
    m_peak[ieta][iphi] = energy;
    m_chi2[ieta][iphi] = chi2;
    h_waveformchi2->Fill(m_peak[ieta][iphi], m_chi2[ieta][iphi]);
    if (m_chi2[ieta][iphi] > 10000)
    {
      m_peak[ieta][iphi] = 0;
    }
    h_time_energy->Fill(time, energy);
  }

  // Apply cut
  for (int ieta = 0; ieta < n_etabin; ++ieta)
  {
    for (int iphi = 0; iphi < n_phibin; ++iphi)
    {
      if (m_peak[ieta][iphi] < tower_threshold)
      {
        continue;  // tower cut
      }
      int up = iphi + 1;
      int down = iphi - 1;
      if (up > 63)
      {
        up -= 64;
      }
      if (down < 0)
      {
        down += 64;
      }
      if (m_peak[ieta][up] < vert_threshold || m_peak[ieta][down] < vert_threshold)
      {
        continue;
      }
      if (ieta != 0 && (m_peak[ieta - 1][up] > veto_threshold || m_peak[ieta - 1][iphi] > veto_threshold || m_peak[ieta - 1][down] > veto_threshold))
      {
        continue;  // left veto cut
      }
      if (ieta != 23 && (m_peak[ieta + 1][up] > veto_threshold || m_peak[ieta + 1][iphi] > veto_threshold || m_peak[ieta + 1][down] > veto_threshold))
      {
        continue;  // right veto cut
      }
      std::cout << "ieta: " << ieta << " iphi: " << iphi << " energy: " << m_peak[ieta][iphi] << " chi2: " << m_chi2[ieta][iphi] << std::endl;
      h_channel_hist[ieta][iphi]->Fill(m_peak[ieta][iphi]);
      h_mip->Fill(m_peak[ieta][iphi]);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int HCalCosmics::ResetEvent(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int HCalCosmics::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "HCalCosmics::End" << std::endl;
  outfile->cd();
  for (auto &ieta : h_channel_hist)
  {
    for (auto &iphi : ieta)
    {
      iphi->Write();
      delete iphi;
    }
  }
  h_mip->Write();
  h_waveformchi2->Write();
  h_time_energy->Write();

  outfile->Close();
  delete outfile;
  return 0;
}

double HCalCosmics::gamma_function(double *x, double *par)
{
  double peak = par[0];
  double shift = par[1];
  double scale = par[2];
  double N = par[3];

  if (scale == 0)
  {
    return 0;
  }

  double arg_para = (x[0] - shift) / scale;
  if (arg_para < 0)
  {
    arg_para = 0;
  }
  double peak_para = (peak - shift) / scale;
  double numerator = N * pow(arg_para, peak_para) * TMath::Exp(-arg_para);
  double denominator = ROOT::Math::tgamma(peak_para + 1) * scale;

  if (denominator == 0)
  {
    return 1e8;
  }
  double val = numerator / denominator;
  if (isnan(val))
  {
    return 0;
  }
  return val;
}

TF1 *HCalCosmics::fitHist(TH1F *h)
{
  TF1 *f_gaus = new TF1("f_gaus", "gaus", 0, 10000);
  h->Fit(f_gaus, "QN", "", 0, 10000);

  TF1 *f_gamma = new TF1("f_gamma", gamma_function, 0, 5000, 4);
  f_gamma->SetParName(0, "Peak(ADC)");
  f_gamma->SetParName(1, "Shift");
  f_gamma->SetParName(2, "Scale");
  f_gamma->SetParName(3, "N");

  f_gamma->SetParLimits(0, 1000, 4000);
  f_gamma->SetParLimits(1, 500, 2000);
  f_gamma->SetParLimits(2, 200, 1000);
  f_gamma->SetParLimits(3, 0, 1e9);

  f_gamma->SetParameter(0, f_gaus->GetParameter(1));
  f_gamma->SetParameter(1, 1300);
  f_gamma->SetParameter(2, 1300);
  f_gamma->SetParameter(3, 1e6);

  h->Fit(f_gamma, "RQN", "", 1000, 5000);

  return f_gamma;
}

void HCalCosmics::fitChannels(const std::string &infile, const std::string &outfilename2)
{
  TFile *fin = new TFile(infile.c_str(), "READ");
  if (!fin)
  {
    std::cout << "file " << infile << "   not found";
    return;
  }
  for (int ieta = 0; ieta < n_etabin; ++ieta)
  {
    for (int iphi = 0; iphi < n_phibin; ++iphi)
    {
      std::string channel_histname = "h_channel_" + std::to_string(ieta) + "_" + std::to_string(iphi);
      h_channel_hist[ieta][iphi] = (TH1F *) fin->Get(channel_histname.c_str());
    }
  }

  if (!h_channel_hist[0][0])
  {
    std::cout << "no hists in " << infile << std::endl;
    return;
  }

  TFile *outfileFit = new TFile(outfilename2.c_str(), "recreate");

  TH2F *h2_peak = new TH2F("h2_peak", "", n_etabin, 0, n_etabin, n_phibin, 0, n_phibin);

  TH1F *h_allTow = (TH1F *) h_channel_hist[0][0]->Clone("h_allTow");
  h_allTow->Reset();

  for (int ieta = 0; ieta < n_etabin; ++ieta)
  {
    for (int iphi = 0; iphi < n_phibin; ++iphi)
    {
      TF1 *res = fitHist(h_channel_hist[ieta][iphi]);
      h2_peak->SetBinContent(ieta + 1, iphi + 1, res->GetParameter(0));
      h_allTow->Add(h_channel_hist[ieta][iphi]);
    }
  }

  TF1 *res = fitHist(h_allTow);
  res->Write("f_allTow");

  h_allTow->Write();
  h2_peak->Write();

  outfileFit->Close();
  fin->Close();
}
