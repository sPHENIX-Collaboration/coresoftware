#include "CentralityCalibrator.h"

#include <cdbobjects/CDBTTree.h>

#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TParameter.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

const std::array<std::string, CentralityCalibrator::NCLASSES> CentralityCalibrator::class_names{"main", "m5", "p5"};
const std::array<std::string, CentralityCalibrator::NCLASSES> CentralityCalibrator::class_suffixes{"", "_m5", "_p5"};

namespace
{
  // vertex bins of the vertex scales: 5 cm for |z| <= 30, 10 cm up to 60, then 30 and 60 cm; the scales are
  // normalised to the 0-5 cm bin
  constexpr int NVERTEXBINS = 22;
  constexpr int REFERENCE_VERTEXBIN = NVERTEXBINS / 2;
  constexpr std::array<double, NVERTEXBINS + 1> vertex_edges{-150, -90, -60, -50, -40, -30, -25, -20, -15, -10, -5, 0,
                                                             5, 10, 15, 20, 25, 30, 40, 50, 60, 90, 150};
  // charge sum histograms
  constexpr int NCHARGEBINS = 2500;
  constexpr double CHARGE_LOW = -0.5;
  constexpr double CHARGE_HIGH = 2499.5;
  // number of rows of the Centrality payload (divisions beyond the calibrated ones repeat the last one)
  constexpr int NPAYLOADDIVS = 100;

  // the events of one class: charge sums of the hit tubes, vertex and ZDC sums
  struct Event
  {
    float ssum;
    float nsum;
    float vertex;
    float zdc_s;
    float zdc_n;
  };

  // keeps the histograms made here out of whatever file is open
  class NoAddDirectory
  {
   public:
    NoAddDirectory()
      : m_status(TH1::AddDirectoryStatus())
    {
      TH1::AddDirectory(false);
    }
    ~NoAddDirectory() { TH1::AddDirectory(m_status); }
    NoAddDirectory(const NoAddDirectory &) = delete;
    NoAddDirectory &operator=(const NoAddDirectory &) = delete;

   private:
    bool m_status;
  };

  // three iterations of a Gaussian fit over mean +- 1.2 sigma around seed: the core of the peak only, so the valleys
  // towards the neighbouring bunches do not pull the mean
  bool coreFit(TH1D &h, double seed, double width0, double &mean, double &err, double &sigma, double &chi2ndf)
  {
    mean = seed;
    sigma = width0;
    TF1 f("fcore", "gaus", seed - (1.2 * width0), seed + (1.2 * width0));
    for (int it = 0; it < 3; it++)
    {
      f.SetRange(mean - (1.2 * sigma), mean + (1.2 * sigma));
      f.SetParameters(h.GetBinContent(h.FindBin(mean)), mean, sigma);
      if (h.Fit(&f, "QRN0L") != 0)
      {
        return false;
      }
      mean = f.GetParameter(1);
      sigma = std::fabs(f.GetParameter(2));
      if (!(sigma > 0.05 && sigma < 3.))
      {
        return false;
      }
    }
    err = f.GetParError(1);
    chi2ndf = f.GetNDF() > 0 ? f.GetChisquare() / f.GetNDF() : -1;
    return true;
  }

  // negative binomial probability of n hits for mean mu and shape k
  double nbd(double n, double mu, double k)
  {
    if (n + k > 100.0)
    {
      // logarithms for large arguments
      const double lngamma = TMath::LnGamma(n + k) - TMath::LnGamma(n + 1.) - TMath::LnGamma(k);
      const double lnpow = (n * TMath::Log(mu / k)) - ((n + k) * TMath::Log(1.0 + (mu / k)));
      return TMath::Exp(lngamma + lnpow);
    }
    return TMath::Gamma(n + k) / (TMath::Gamma(n + 1.) * TMath::Gamma(k)) *
           TMath::Exp((n * TMath::Log(mu / k)) - ((n + k) * TMath::Log(1.0 + (mu / k))));
  }

  struct NbdGlauberFit
  {
    double mu{0};
    double k{0};
    double chi2{0};
    double trigeff{0};
    double trigeff_err{0};
  };

  // NBD+Glauber fit of a min-bias charge sum: each participant contributes an NBD(mu, k), summed over the Glauber
  // Npart distribution, fitted over the fully efficient bins 700-2000. The trigger efficiency is the data integral over
  // that of the fitted model normalised to the data above bin 700. hsim gets the model.
  NbdGlauberFit fitNbdGlauber(const TH1D &hdata, const TH1D &hnpart, TH1D &hsim)
  {
    constexpr int lowfit = 700;
    constexpr int highfit = 2000;
    const double cut = hdata.Integral(lowfit, highfit);
    const double all = hdata.Integral();
    const double after = cut / all;
    TH1D hscaled(hdata);
    hscaled.SetName("h_charge_sum_nbdfit");
    hscaled.Scale(hnpart.Integral() * after / all);

    const double alpha = 1.;
    TF1 fnbd(
        "nbd", [&hnpart](const double *x, const double *par)
        {
          double result = 0;
          for (int ib = 2; ib <= hnpart.GetNbinsX(); ib++)
          {
            const int npart = static_cast<int>(hnpart.GetBinCenter(ib));
            const double weight = hnpart.GetBinContent(ib);
            if (weight <= 0)
            {
              continue;
            }
            result += nbd(x[0], par[0] * std::pow(npart, par[2]), par[1] * std::pow(npart, par[2])) * weight;
          }
          return par[3] * result; },
        0, 2500, 4);
    fnbd.SetParNames("mu", "k", "alpha", "N");

    // the chi2 surface has local minima: the fit starts from a grid of (mu, k) and keeps the best one
    NbdGlauberFit fit;
    fit.chi2 = 1e300;
    for (double mu0 : {4.7, 1.0, 2.0, 4.0, 7.0, 10.0})
    {
      for (double k0 : {1.16, 0.5, 2.0})
      {
        fnbd.SetParameters(mu0, k0, alpha, 1.0);
        fnbd.SetParLimits(0, 1, 10);
        fnbd.SetParLimits(1, 0.1, 3);
        fnbd.FixParameter(2, alpha);
        fnbd.SetParLimits(3, 0.1, 1.5);
        hscaled.Fit(&fnbd, "NDOQR", "", lowfit, highfit);
        if (fnbd.GetChisquare() < fit.chi2)
        {
          fit.chi2 = fnbd.GetChisquare();
          fit.mu = fnbd.GetParameter(0);
          fit.k = fnbd.GetParameter(1);
        }
      }
    }

    hsim.Reset();
    for (int ib = 2; ib <= hnpart.GetNbinsX(); ib++)
    {
      const int npart = static_cast<int>(hnpart.GetBinCenter(ib));
      const double weight = hnpart.GetBinContent(ib);
      if (weight <= 0)
      {
        continue;
      }
      for (int ihit = 0; ihit < (20 * npart * static_cast<int>(fit.mu)); ihit++)
      {
        hsim.Fill(static_cast<double>(ihit), nbd(ihit, fit.mu * std::pow(npart, alpha), fit.k * std::pow(npart, alpha)) * weight);
      }
    }
    const int nbins = hdata.GetNbinsX();
    hsim.Scale(hdata.Integral(lowfit, nbins) / hsim.Integral(lowfit, nbins));
    double err_data;
    double err_sim;
    const double integral_data = hdata.IntegralAndError(1, nbins, err_data);
    const double integral_sim = hsim.IntegralAndError(1, nbins, err_sim);
    fit.trigeff = integral_data / integral_sim;
    fit.trigeff_err = fit.trigeff * std::sqrt(TMath::Power(err_data / hdata.Integral(1, nbins), 2) + TMath::Power(err_sim / hsim.Integral(1, nbins), 2));
    return fit;
  }

  // the per-class tables a payload is made of, read from a Calibrate output file
  struct VertexBin
  {
    float low{0};
    float high{0};
    std::array<float, CentralityCalibrator::NCLASSES> scale{};
  };

  bool readDivisions(TFile &f, std::array<std::vector<float>, CentralityCalibrator::NCLASSES> &divs)
  {
    for (int c = 0; c < CentralityCalibrator::NCLASSES; c++)
    {
      TNtuple *tn = dynamic_cast<TNtuple *>(f.Get((CentralityCalibrator::class_names[c] + "/tn_centrality").c_str()));
      if (!tn || tn->GetEntries() < 1)
      {
        return false;
      }
      float low;
      tn->SetBranchAddress("low", &low);
      divs[c].clear();
      for (Long64_t i = 0; i < tn->GetEntries(); i++)
      {
        tn->GetEntry(i);
        divs[c].push_back(low);
      }
    }
    return true;
  }

  bool readScales(TFile &f, std::array<double, CentralityCalibrator::NCLASSES> &scales)
  {
    for (int c = 0; c < CentralityCalibrator::NCLASSES; c++)
    {
      TNtuple *tn = dynamic_cast<TNtuple *>(f.Get((CentralityCalibrator::class_names[c] + "/tn_scale").c_str()));
      if (!tn || tn->GetEntries() < 1)
      {
        return false;
      }
      float scale;
      tn->SetBranchAddress("scale", &scale);
      tn->GetEntry(0);
      scales[c] = scale;
    }
    return true;
  }

  // non-finite scales (empty vertex bins outside the vertex cut) become 0, which is also what the reco uses outside
  // every vertex bin
  bool readVertexScales(TFile &f, std::vector<VertexBin> &bins)
  {
    bins.clear();
    for (int c = 0; c < CentralityCalibrator::NCLASSES; c++)
    {
      TNtuple *tn = dynamic_cast<TNtuple *>(f.Get((CentralityCalibrator::class_names[c] + "/tn_vertexscale").c_str()));
      if (!tn)
      {
        return false;
      }
      float low;
      float high;
      float scale;
      tn->SetBranchAddress("low_vertex", &low);
      tn->SetBranchAddress("high_vertex", &high);
      tn->SetBranchAddress("scale", &scale);
      if (c == 0)
      {
        bins.resize(tn->GetEntries());
      }
      if (static_cast<Long64_t>(bins.size()) != tn->GetEntries())
      {
        return false;
      }
      for (Long64_t i = 0; i < tn->GetEntries(); i++)
      {
        tn->GetEntry(i);
        if (c == 0)
        {
          bins[i].low = low;
          bins[i].high = high;
        }
        else if (low != bins[i].low || high != bins[i].high)
        {
          return false;  // the classes share the vertex bins
        }
        bins[i].scale[c] = std::isfinite(scale) ? scale : 0.F;
      }
    }
    return true;
  }

  template <class T>
  bool readParameter(TFile &f, const char *name, T &value)
  {
    TParameter<T> *p = dynamic_cast<TParameter<T> *>(f.Get(name));
    if (!p)
    {
      return false;
    }
    value = p->GetVal();
    return true;
  }
}  // namespace

int CentralityCalibrator::ClassOf(double t0) const
{
  // windows centred on each bunch, edges halfway to the neighbouring bunches; the inner edges belong to the main
  // bunch; NaN is in no class
  const double hm = 0.5 * m_spacing_m;
  const double hp = 0.5 * m_spacing_p;
  if (t0 >= -hm && t0 <= hp)
  {
    return 0;
  }
  if (t0 >= -3 * hm && t0 < -hm)
  {
    return 1;
  }
  if (t0 > hp && t0 <= 3 * hp)
  {
    return 2;
  }
  return -1;
}

bool CentralityCalibrator::ReadT0Offset(const std::string &t0table, int run, double &offset, int &status)
{
  std::ifstream in(t0table);
  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line))
  {
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> v;
    while (std::getline(ss, token, ','))
    {
      v.push_back(token);
    }
    if (v.size() < 11 || std::atoi(v[0].c_str()) != run)
    {
      continue;
    }
    offset = std::atof(v[2].c_str());
    status = std::atoi(v[10].c_str());
    return true;
  }
  return false;
}

int CentralityCalibrator::FitT0Offset(int run, const std::string &treefile, const std::string &outdir) const
{
  TFile *fin = TFile::Open(treefile.c_str(), "READ");
  TTree *tree = (fin && !fin->IsZombie()) ? dynamic_cast<TTree *>(fin->Get("T")) : nullptr;
  if (!tree)
  {
    std::cout << "CentralityCalibrator::FitT0Offset: no tree T in " << treefile << std::endl;
    delete fin;
    return 1;
  }
  ULong64_t gl1_scaled;
  float zdc_sum[2];
  float tube_charge[128];
  float tube_time[128];
  float vertex;
  float t0;
  tree->SetBranchStatus("*", false);
  for (const char *branch : {"gl1_scaled", "zdc_sum", "mbd_charge", "mbd_time", "mbd_vertex", "mbd_time_zero"})
  {
    tree->SetBranchStatus(branch, true);
  }
  tree->SetBranchAddress("gl1_scaled", &gl1_scaled);
  tree->SetBranchAddress("zdc_sum", zdc_sum);
  tree->SetBranchAddress("mbd_charge", tube_charge);
  tree->SetBranchAddress("mbd_time", tube_time);
  tree->SetBranchAddress("mbd_vertex", &vertex);
  tree->SetBranchAddress("mbd_time_zero", &t0);

  NoAddDirectory noadd;
  // 0.05 ns bins; the main peak sigma is ~0.8 ns
  TH1D h("h_t0", ";mbd_time_zero [ns];events", 480, -12, 12);
  long nsel = 0;
  for (Long64_t i = 0; i < tree->GetEntries() && nsel < m_t0fit_max_events; i++)
  {
    tree->GetEntry(i);
    if (!((gl1_scaled >> m_trigger_bit) & 0x1U))
    {
      continue;
    }
    int nhit_s = 0;
    int nhit_n = 0;
    for (int ich = 0; ich < 64; ich++)
    {
      if (std::fabs(tube_time[ich]) < m_hit_time_cut && tube_charge[ich] > m_hit_charge_cut)
      {
        nhit_s++;
      }
      if (std::fabs(tube_time[ich + 64]) < m_hit_time_cut && tube_charge[ich + 64] > m_hit_charge_cut)
      {
        nhit_n++;
      }
    }
    if (nhit_s < m_nhit_cut || nhit_n < m_nhit_cut || zdc_sum[0] <= m_zdc_cut || zdc_sum[1] <= m_zdc_cut ||
        std::fabs(vertex) > m_t0fit_vertex_cut)
    {
      continue;
    }
    nsel++;
    h.Fill(t0);
  }
  delete fin;

  // status: 0 ok; 1 fewer than 1000 events (offset 0); 2 main-bunch fit failed; 3 |offset| > 3 ns (used, but check
  // the run); 4 ambiguous main bunch: the tallest bin over the full range is outside the seed window, i.e. a side
  // bunch is taller than the bunch at the nominal time (the offset is for information only)
  int status = 0;
  double offset = 0;
  double offset_err = 0;
  double sigma = 0;
  double chi2ndf = -1;
  std::array<double, 2> side_rel{-99, -99};  // -5 and +5 ns peaks relative to the main one, -99 = not found
  std::array<double, 2> side_err{0, 0};
  if (nsel < 1000)
  {
    status = 1;
  }
  else
  {
    h.GetXaxis()->SetRangeUser(-m_t0fit_seed_window, m_t0fit_seed_window);
    const double seed = h.GetBinCenter(h.GetMaximumBin());
    h.GetXaxis()->SetRange();
    const bool ambiguous = std::fabs(h.GetBinCenter(h.GetMaximumBin())) >= m_t0fit_seed_window;
    if (!coreFit(h, seed, 0.8, offset, offset_err, sigma, chi2ndf))
    {
      status = 2;
      offset = 0;
    }
    else
    {
      if (std::fabs(offset) > 3.)
      {
        status = 3;
      }
      if (ambiguous)
      {
        status = 4;
      }
      const std::array<double, 2> side_seed{offset - 5.1, offset + 5.1};
      for (int s = 0; s < 2; s++)
      {
        // seeded at the tallest bin within +-1.5 ns of the expected position
        h.GetXaxis()->SetRangeUser(side_seed[s] - 1.5, side_seed[s] + 1.5);
        const double sseed = h.GetBinCenter(h.GetMaximumBin());
        h.GetXaxis()->SetRange();
        double mean;
        double err;
        double sig;
        double chi2;
        if (h.Integral(h.FindBin(sseed - 1.), h.FindBin(sseed + 1.)) > 500 && coreFit(h, sseed, 0.8, mean, err, sig, chi2) &&
            std::fabs(mean - side_seed[s]) < 1.5)
        {
          side_rel[s] = mean - offset;
          side_err[s] = err;
        }
      }
    }
  }

  std::filesystem::create_directories(outdir);
  const std::string stem = outdir + "/t0offset_" + std::to_string(run);
  FILE *fout = fopen((stem + ".csv").c_str(), "w");
  if (!fout)
  {
    std::cout << "CentralityCalibrator::FitT0Offset: cannot write " << stem << ".csv" << std::endl;
    return 1;
  }
  fprintf(fout, "%d,%ld,%.4f,%.4f,%.4f,%.3f,%.4f,%.4f,%.4f,%.4f,%d\n", run, nsel, offset, offset_err, sigma, chi2ndf,
          side_rel[0], side_err[0], side_rel[1], side_err[1], status);
  fclose(fout);
  std::cout << "run " << run << ": " << nsel << " events, t0 offset " << offset << " +- " << offset_err << " ns, sigma "
            << sigma << ", -5/+5 ns peaks at " << side_rel[0] << " / " << side_rel[1] << " ns, status " << status << std::endl;

  TFile fhist((stem + ".root").c_str(), "RECREATE");
  h.Write();
  fhist.Close();
  return 0;
}

int CentralityCalibrator::Calibrate(int run, const std::string &treefile, const std::string &t0table, const std::string &reffile,
                                    const std::string &outfile) const
{
  double t0offset = 0;
  int t0status = 9;
  const bool isref = reffile.empty();
  if (!ReadT0Offset(t0table, run, t0offset, t0status) || (t0status != 0 && t0status != 3))
  {
    // a run without a usable main-bunch offset is not calibrated: its payloads get the reference run's tables
    // (WritePayloads), so that is not a failure, unless the run is a reference run
    std::cout << "CentralityCalibrator::Calibrate: run " << run << " has no usable t0 offset in " << t0table << " (status "
              << t0status << "): " << (isref ? "it cannot be a reference run" : "not calibrated, its payloads use the reference run's tables")
              << std::endl;
    return isref ? 1 : 0;
  }

  NoAddDirectory noadd;
  // the reference run's mean charge sums, unbalanced and balanced, per class
  int refrun = run;
  std::array<std::array<float, 2>, NCLASSES> refmean{};
  if (!isref)
  {
    TFile fref(reffile.c_str(), "READ");
    bool ok = !fref.IsZombie() && readParameter(fref, "run", refrun);
    for (int c = 0; c < NCLASSES && ok; c++)
    {
      TH1D *h = dynamic_cast<TH1D *>(fref.Get((class_names[c] + "/h_charge_sum").c_str()));
      TH1D *hbal = dynamic_cast<TH1D *>(fref.Get((class_names[c] + "/h_charge_sum_balanced").c_str()));
      ok = h && hbal;
      if (ok)
      {
        refmean[c] = {static_cast<float>(h->GetMean()), static_cast<float>(hbal->GetMean())};
      }
    }
    if (!ok)
    {
      std::cout << "CentralityCalibrator::Calibrate: " << reffile << " is not a reference calibration" << std::endl;
      return 1;
    }
  }

  // reference run: the Glauber Npart distribution of the NBD+Glauber fit
  TH1D *hnpart = nullptr;
  if (isref)
  {
    TFile fglauber(m_glauber_file.c_str(), "READ");
    TH1D *h = fglauber.IsZombie() ? nullptr : dynamic_cast<TH1D *>(fglauber.Get("hNpart"));
    if (!h)
    {
      std::cout << "CentralityCalibrator::Calibrate: the reference run needs a Glauber file with hNpart (SetGlauberFile), got \""
                << m_glauber_file << "\"" << std::endl;
      return 1;
    }
    hnpart = dynamic_cast<TH1D *>(h->Clone("hNpart"));
  }
  std::unique_ptr<TH1D> hnpart_owner(hnpart);

  TFile *fin = TFile::Open(treefile.c_str(), "READ");
  TTree *tree = (fin && !fin->IsZombie()) ? dynamic_cast<TTree *>(fin->Get("T")) : nullptr;
  if (!tree)
  {
    std::cout << "CentralityCalibrator::Calibrate: no tree T in " << treefile << std::endl;
    delete fin;
    return 1;
  }
  ULong64_t gl1_scaled;
  float zdc_sum[2];
  float tube_charge[128];
  float tube_time[128];
  float vertex;
  float t0;
  tree->SetBranchStatus("*", false);
  for (const char *branch : {"gl1_scaled", "zdc_sum", "mbd_charge", "mbd_time", "mbd_vertex", "mbd_time_zero"})
  {
    tree->SetBranchStatus(branch, true);
  }
  tree->SetBranchAddress("gl1_scaled", &gl1_scaled);
  tree->SetBranchAddress("zdc_sum", zdc_sum);
  tree->SetBranchAddress("mbd_charge", tube_charge);
  tree->SetBranchAddress("mbd_time", tube_time);
  tree->SetBranchAddress("mbd_vertex", &vertex);
  tree->SetBranchAddress("mbd_time_zero", &t0);

  // one pass over the run: triggered events with enough hit tubes, sorted into their class
  std::array<std::vector<Event>, NCLASSES> events;
  long ntrig = 0;
  for (Long64_t i = 0; i < tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    t0 -= t0offset;
    if (!((gl1_scaled >> m_trigger_bit) & 0x1U))
    {
      continue;
    }
    ntrig++;
    int nhit_s = 0;
    int nhit_n = 0;
    float ssum = 0;
    float nsum = 0;
    for (int ich = 0; ich < 64; ich++)
    {
      if (std::fabs(tube_time[ich] - t0offset) < m_hit_time_cut && tube_charge[ich] > m_hit_charge_cut)
      {
        nhit_s++;
        ssum += tube_charge[ich];
      }
      if (std::fabs(tube_time[ich + 64] - t0offset) < m_hit_time_cut && tube_charge[ich + 64] > m_hit_charge_cut)
      {
        nhit_n++;
        nsum += tube_charge[ich + 64];
      }
    }
    if (nhit_s < m_nhit_cut || nhit_n < m_nhit_cut)
    {
      continue;
    }
    const int c = ClassOf(t0);
    if (c >= 0)
    {
      events[c].push_back({ssum, nsum, vertex, zdc_sum[0], zdc_sum[1]});
    }
  }
  const Long64_t nentries = tree->GetEntries();
  delete fin;
  std::cout << "CentralityCalibrator::Calibrate: run " << run << " (reference " << refrun << "), t0 offset " << t0offset << " ns: "
            << nentries << " events, " << ntrig << " with trigger bit " << m_trigger_bit << std::endl;

  std::filesystem::path outpath(outfile);
  if (outpath.has_parent_path())
  {
    std::filesystem::create_directories(outpath.parent_path());
  }
  TFile fout(outfile.c_str(), "RECREATE");
  TParameter<int>("run", run).Write();
  TParameter<int>("refrun", refrun).Write();
  TParameter<int>("trigger_bit", m_trigger_bit).Write();
  TParameter<double>("t0offset", t0offset).Write();
  TParameter<int>("t0status", t0status).Write();
  TParameter<double>("spacing_m", m_spacing_m).Write();
  TParameter<double>("spacing_p", m_spacing_p).Write();

  for (int c = 0; c < NCLASSES; c++)
  {
    // min-bias events of the class: ZDC coincidence and vertex cut
    TH2D h2("h2_charge_sum_vs_vertex", ";z_{vtx} [cm];charge sum", NVERTEXBINS, vertex_edges.data(), NCHARGEBINS, CHARGE_LOW, CHARGE_HIGH);
    TH1D h("h_charge_sum", ";charge sum;events", NCHARGEBINS, CHARGE_LOW, CHARGE_HIGH);
    TH1D hbal("h_charge_sum_balanced", ";vertex-balanced charge sum;events", NCHARGEBINS, CHARGE_LOW, CHARGE_HIGH);
    std::vector<Event> minbias;
    for (const Event &ev : events[c])
    {
      if (ev.zdc_s <= m_zdc_cut || ev.zdc_n <= m_zdc_cut || std::fabs(ev.vertex) > m_vertex_cut)
      {
        continue;
      }
      h2.Fill(ev.vertex, ev.ssum + ev.nsum);
      h.Fill(ev.ssum + ev.nsum);
      minbias.push_back(ev);
    }

    // vertex scales: mean charge sum of the 0-5 cm bin over that of each bin
    std::array<double, NVERTEXBINS> vertex_scale{};
    for (int ib = 0; ib < NVERTEXBINS; ib++)
    {
      TH1D *proj = h2.ProjectionY(Form("h_charge_sum_vertexbin_%d", ib), ib + 1, ib + 1);
      vertex_scale[ib] = proj->GetMean();
      delete proj;
    }
    const double reference_mean = vertex_scale[REFERENCE_VERTEXBIN];
    for (double &scale : vertex_scale)
    {
      scale = reference_mean / scale;
    }
    for (const Event &ev : minbias)
    {
      double scale = 0;
      for (int iv = 0; iv < NVERTEXBINS; iv++)
      {
        if (ev.vertex < vertex_edges[iv + 1])
        {
          scale = vertex_scale[iv];
          break;
        }
      }
      hbal.Fill((ev.ssum * scale) + (ev.nsum * scale));
    }

    // run scale: the reference run's balanced mean over this run's
    const float run_scale = isref ? 1.F : static_cast<float>(refmean[c][1] / hbal.GetMean());
    std::cout << "  class " << class_names[c] << ": " << events[c].size() << " events in the window, " << minbias.size()
              << " min-bias; mean charge sum " << h.GetMean() << ", balanced " << hbal.GetMean() << ", run scale " << run_scale
              << std::endl;

    TDirectory *dir = fout.mkdir(class_names[c].c_str());
    dir->cd();
    TParameter<Long64_t>("nevents_window", static_cast<Long64_t>(events[c].size())).Write();
    TParameter<Long64_t>("nevents_minbias", static_cast<Long64_t>(minbias.size())).Write();
    h2.Write();
    h.Write();
    hbal.Write();
    TNtuple tn_vertexscale("tn_vertexscale", "vertex scales", "bin:low_vertex:high_vertex:scale");
    for (int iv = 0; iv < NVERTEXBINS; iv++)
    {
      tn_vertexscale.Fill(iv, vertex_edges[iv], vertex_edges[iv + 1], vertex_scale[iv]);
    }
    tn_vertexscale.Write();
    TNtuple tn_scale("tn_scale", "run scale", "run:scale");
    tn_scale.Fill(run, run_scale);
    tn_scale.Write();

    if (isref)
    {
      // NBD+Glauber fits of the charge sum and of the balanced one: the trigger efficiency
      for (TH1D *hfit : {&h, &hbal})
      {
        const bool balanced = (hfit == &hbal);
        TH1D hsim(*hfit);
        hsim.SetName(balanced ? "hSimMBD_balanced" : "hSimMBD");
        const NbdGlauberFit fit = fitNbdGlauber(*hfit, *hnpart, hsim);
        const std::string sfx = balanced ? "_balanced" : "";
        std::cout << "  class " << class_names[c] << (balanced ? ", balanced" : "") << ": NBD+Glauber mu " << fit.mu << " k " << fit.k
                  << " chi2 " << fit.chi2 << ", trigger efficiency " << fit.trigeff << " +- " << fit.trigeff_err << std::endl;
        hsim.Write();
        TParameter<double>(("nbd_mu" + sfx).c_str(), fit.mu).Write();
        TParameter<double>(("nbd_k" + sfx).c_str(), fit.k).Write();
        TParameter<double>(("nbd_chi2" + sfx).c_str(), fit.chi2).Write();
        TParameter<double>(("trigeff" + sfx).c_str(), fit.trigeff).Write();
        TParameter<double>(("trigeff_err" + sfx).c_str(), fit.trigeff_err).Write();
      }

      // divisions: the balanced charge sums (vertex scales as stored) sorted from the most central event, in
      // m_ndivisions equal-population bins; bin i spans [low, high]
      std::vector<float> sums;
      sums.reserve(minbias.size());
      for (const Event &ev : minbias)
      {
        float scale = 0;
        bool found = false;
        for (int iv = 0; iv < NVERTEXBINS && !found; iv++)
        {
          const auto low = static_cast<float>(vertex_edges[iv]);
          const auto high = static_cast<float>(vertex_edges[iv + 1]);
          if (ev.vertex < high && ev.vertex >= low)
          {
            scale = static_cast<float>(vertex_scale[iv]);
            found = true;
          }
        }
        const double total = (ev.ssum + ev.nsum) * static_cast<double>(scale);
        if (!found || total > m_max_charge_sum)
        {
          continue;
        }
        sums.push_back(static_cast<float>(total));
      }
      std::sort(sums.begin(), sums.end(), [](float a, float b) { return a > b; });
      const int size = static_cast<int>(sums.size());
      if (size < m_ndivisions)
      {
        std::cout << "CentralityCalibrator::Calibrate: class " << class_names[c] << " has " << size << " events for "
                  << m_ndivisions << " divisions" << std::endl;
        return 1;
      }
      TNtuple tn_centrality("tn_centrality", "centrality divisions", "bin:low:high");
      for (int i = 0; i < m_ndivisions; i++)
      {
        // float arithmetic as in the original calibration, which the existing payloads come from
        int ibin = static_cast<int>(std::floor(static_cast<float>(i) * static_cast<float>(size) / static_cast<float>(m_ndivisions)));
        int ibin1 = static_cast<int>(std::floor(static_cast<float>(i + 1) * static_cast<float>(size) / static_cast<float>(m_ndivisions)));
        if (i == m_ndivisions - 1)
        {
          ibin1 = size - 1;
        }
        ibin1 = std::min(ibin1, size - 1);
        tn_centrality.Fill(i, sums[ibin1], sums[ibin]);
      }
      tn_centrality.Write();
      TParameter<int>("ndivisions", m_ndivisions).Write();
    }
    fout.cd();
  }
  fout.Close();
  std::cout << "CentralityCalibrator::Calibrate: wrote " << outfile << std::endl;
  return 0;
}

int CentralityCalibrator::WritePayloads(int run, const std::string &calibfile, const std::string &reffile, const std::string &t0table,
                                        const std::string &outdir)
{
  // reference run: divisions, vertex scales and the bunch spacings
  std::array<std::vector<float>, NCLASSES> divs;
  std::vector<VertexBin> refvertex;
  int refrun = 0;
  double spacing_m = 0;
  double spacing_p = 0;
  {
    TFile fref(reffile.c_str(), "READ");
    if (fref.IsZombie() || !readParameter(fref, "run", refrun) || !readParameter(fref, "spacing_m", spacing_m) ||
        !readParameter(fref, "spacing_p", spacing_p) || !readDivisions(fref, divs) || !readVertexScales(fref, refvertex))
    {
      std::cout << "CentralityCalibrator::WritePayloads: " << reffile << " is not a reference calibration" << std::endl;
      return -1;
    }
  }

  // the run's own calibration, if it has a usable one: t0 status 0 or 3, finite positive scales
  double t0offset = 0;
  int t0status = 9;
  if (!ReadT0Offset(t0table, run, t0offset, t0status))
  {
    t0status = 9;
  }
  std::array<double, NCLASSES> scales{};
  std::vector<VertexBin> vertex;
  bool own = false;
  if (!calibfile.empty() && (t0status == 0 || t0status == 3) && std::filesystem::exists(calibfile))
  {
    TFile fcal(calibfile.c_str(), "READ");
    own = !fcal.IsZombie() && readScales(fcal, scales) && readVertexScales(fcal, vertex);
    for (double scale : scales)
    {
      own = own && std::isfinite(scale) && scale > 0;
    }
  }
  PayloadSource source = OWN;
  if (run == refrun)
  {
    source = REFERENCE_RUN;
  }
  else if (!own)
  {
    source = REFERENCE_DEFAULT;
  }
  if (source != OWN)
  {
    scales.fill(1.);
    vertex = refvertex;
  }
  // the fitted offset where there is one (status 4: fitted but ambiguous), else 0
  const double offset = (t0status == 0 || t0status == 3 || t0status == 4) ? t0offset : 0.;

  for (const char *sub : {"divs", "scales", "vertexscales"})
  {
    std::filesystem::create_directories(outdir + "/" + sub);
  }
  const std::string srun = std::to_string(run);
  {
    CDBTTree cdbttree(outdir + "/divs/cdb_centrality_" + srun + ".root");
    for (int c = 0; c < NCLASSES; c++)
    {
      const int ndivs = static_cast<int>(divs[c].size());
      for (int i = 0; i < NPAYLOADDIVS; i++)
      {
        cdbttree.SetFloatValue(i, "centralitydiv" + class_suffixes[c], divs[c][std::min(i, ndivs - 1)]);
      }
    }
    cdbttree.Commit();
    cdbttree.WriteCDBTTree();
  }
  {
    CDBTTree cdbttree(outdir + "/scales/cdb_centrality_scale_" + srun + ".root");
    for (int c = 0; c < NCLASSES; c++)
    {
      cdbttree.SetDoubleValue(0, "centralityscale" + class_suffixes[c], scales[c]);
    }
    cdbttree.SetDoubleValue(0, "t0offset", offset);
    cdbttree.SetIntValue(0, "t0status", t0status);
    cdbttree.SetDoubleValue(0, "spacing_m", spacing_m);
    cdbttree.SetDoubleValue(0, "spacing_p", spacing_p);
    cdbttree.Commit();
    cdbttree.WriteCDBTTree();
  }
  {
    CDBTTree cdbttree(outdir + "/vertexscales/cdb_centrality_vertex_scale_" + srun + ".root");
    cdbttree.SetIntValue(0, "nvertexbins", static_cast<int>(vertex.size()));
    for (size_t i = 0; i < vertex.size(); i++)
    {
      const int channel = static_cast<int>(i);
      cdbttree.SetDoubleValue(channel, "low_vertex", vertex[i].low);
      cdbttree.SetDoubleValue(channel, "high_vertex", vertex[i].high);
      for (int c = 0; c < NCLASSES; c++)
      {
        cdbttree.SetDoubleValue(channel, "scale" + class_suffixes[c], vertex[i].scale[c]);
      }
    }
    cdbttree.Commit();
    cdbttree.WriteCDBTTree();
  }
  return source;
}
