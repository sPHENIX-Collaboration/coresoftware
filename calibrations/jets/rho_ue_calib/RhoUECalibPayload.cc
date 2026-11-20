#include "RhoUECalibPayload.h"

#include "RhoUECalibBase.h"
#include "RhoUECalibMaker.h"

#include <cdbobjects/CDBHistos.h>
#include <cdbobjects/CDBTTree.h>

#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TKey.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  double edge_func(const double *xx, const double *par)
  {
    const double dq = xx[0] - par[2];
    return (par[0] + (par[1] * dq)) * 0.5 * std::erfc(dq / (std::sqrt(2.) * par[3]));
  }

  std::unique_ptr<TFile> open_file(const std::string &name)
  {
    std::unique_ptr<TFile> f(TFile::Open(name.c_str(), "READ"));
    if (!f || f->IsZombie())
    {
      std::cout << "RhoUECalibPayload: cannot open " << name << std::endl;
      return nullptr;
    }
    return f;
  }
}  // namespace

RhoUECalibPayload::Knee RhoUECalibPayload::FitKnee(TH1 *hq)
{
  Knee k;
  if (hq->Integral() < 1000)
  {
    return k;
  }
  double prob = 0.995;
  double q995 = 0;
  hq->GetQuantiles(1, &q995, &prob);
  // about 100 bins up to the edge
  std::unique_ptr<TH1> h(static_cast<TH1 *>(hq->Clone("knee_fit")));
  h->SetDirectory(nullptr);
  h->Rebin(std::max(1, static_cast<int>(std::lround(q995 / 100. / hq->GetBinWidth(1)))));
  double plateau = 0;
  int nplat = 0;
  for (int b = h->FindBin(0.6 * q995); b <= h->FindBin(0.8 * q995); ++b)
  {
    plateau += h->GetBinContent(b);
    ++nplat;
  }
  plateau /= std::max(1, nplat);
  const double slope = (h->GetBinContent(h->FindBin(0.85 * q995)) - h->GetBinContent(h->FindBin(0.6 * q995))) / (0.25 * q995);
  TF1 f("knee_edge", edge_func, 0.55 * q995, 1.35 * q995, 4);
  f.SetParameters(plateau, slope, 0.97 * q995, 0.04 * q995);
  f.SetParLimits(2, 0.6 * q995, 1.3 * q995);
  f.SetParLimits(3, 0.002 * q995, 0.5 * q995);
  k.status = h->Fit(&f, "RLQ0N");
  k.knee = f.GetParameter(2);
  k.knee_err = f.GetParError(2);
  k.sigma = f.GetParameter(3);
  k.chi2ndf = f.GetNDF() > 0 ? f.GetChisquare() / f.GetNDF() : -1;
  return k;
}

int RhoUECalibPayload::FitMbdKnees(const std::string &sums_file, const std::string &table, const double t0_window)
{
  auto f = open_file(sums_file);
  if (!f)
  {
    return 1;
  }
  std::ofstream out(table);
  out << "# MBD charge knee and main-bunch t0 centre per run (RhoUECalibPayload::FitMbdKnees, " << sums_file << ")\n"
      << "# run knee knee_err sigma chi2ndf status events t0 t0_rms\n";
  int nfit = 0;
  int nbad = 0;
  for (auto *obj : *f->GetListOfKeys())
  {
    const std::string name = obj->GetName();
    if (name.rfind("q_", 0) != 0)
    {
      continue;
    }
    auto *h = f->Get<TH1>(name.c_str());
    const Knee k = FitKnee(h);
    const bool ok = (k.status == 0 && k.knee > 0);
    // main-bunch t0: the mode, then the mean within +- t0_window, iterated
    double t0 = std::nan("");
    double t0_rms = std::nan("");
    auto *ht0 = f->Get<TH1>(("t0_" + name.substr(2)).c_str());
    if (ht0 && ht0->Integral() > 0)
    {
      t0 = ht0->GetBinCenter(ht0->GetMaximumBin());
      for (int iter = 0; iter < 5; ++iter)
      {
        ht0->GetXaxis()->SetRangeUser(t0 - t0_window, t0 + t0_window);
        t0 = ht0->GetMean();
      }
      t0_rms = ht0->GetRMS();
      ht0->GetXaxis()->SetRange(0, 0);
    }
    out << name.substr(2) << " " << (ok ? k.knee : -1) << " " << k.knee_err << " " << k.sigma << " " << k.chi2ndf << " "
        << k.status << " " << h->Integral() << " " << t0 << " " << t0_rms << "\n";
    ++nfit;
    nbad += ok ? 0 : 1;
  }
  std::cout << "RhoUECalibPayload::FitMbdKnees - " << nfit << " runs, " << nbad << " failed (knee -1) -> " << table << std::endl;
  return 0;
}

int RhoUECalibPayload::MakeKneePayloads(const std::string &table, const std::string &out_dir, const double min_events,
                                        const double max_sigma_dev)
{
  std::ifstream in(table);
  if (!in)
  {
    std::cout << "RhoUECalibPayload::MakeKneePayloads - cannot read " << table << std::endl;
    return 1;
  }
  struct Row
  {
    int run{0};
    double knee{-1};
    double sigma{-1};
    double events{0};
  };
  std::vector<Row> rows;
  std::vector<double> knees;
  std::vector<double> sigmas;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#')
    {
      continue;
    }
    std::istringstream ss(line);
    Row r;
    double knee_err = 0;
    double chi2ndf = 0;
    int status = 0;
    if (ss >> r.run >> r.knee >> knee_err >> r.sigma >> chi2ndf >> status >> r.events)
    {
      rows.push_back(r);
      if (r.knee > 0)
      {
        knees.push_back(r.knee);
        sigmas.push_back(r.sigma);
      }
    }
  }
  if (knees.empty())
  {
    std::cout << "RhoUECalibPayload::MakeKneePayloads - no fitted knee in " << table << std::endl;
    return 1;
  }
  auto median = [](std::vector<double> v)
  {
    std::nth_element(v.begin(), v.begin() + (v.size() / 2), v.end());
    return v[v.size() / 2];
  };
  const double knee_med = median(knees);
  const double sigma_med = median(sigmas);

  std::filesystem::create_directories(out_dir);
  std::ofstream iov(out_dir + "/mbdq_knee_iov.txt");
  iov << "# per-run MBD charge knee payloads (RhoUECalibPayload::MakeKneePayloads, " << table << "), each valid for its run\n"
      << "# run status mbdq_knee file (status 1: the median knee " << knee_med << ")\n";
  int nsub = 0;
  for (const auto &r : rows)
  {
    const bool good = r.knee > 0 && r.events >= min_events && std::abs((r.sigma / sigma_med) - 1) <= max_sigma_dev;
    const double knee = good ? r.knee : knee_med;
    const std::string file = out_dir + "/mbdq_knee_" + std::to_string(r.run) + ".root";
    CDBTTree cdb(file);
    cdb.SetSingleFloatValue("mbdq_knee", static_cast<float>(knee));
    cdb.SetSingleFloatValue("mbdq_knee_fit", static_cast<float>(r.knee));
    cdb.SetSingleFloatValue("mbdq_knee_sigma", static_cast<float>(r.sigma));
    cdb.SetSingleFloatValue("events", static_cast<float>(r.events));
    cdb.SetSingleIntValue("run", r.run);
    cdb.SetSingleIntValue("status", good ? 0 : 1);
    cdb.CommitSingle();
    cdb.WriteCDBTTree();
    iov << r.run << " " << (good ? 0 : 1) << " " << knee << " " << file << "\n";
    if (!good)
    {
      ++nsub;
      std::cout << "RhoUECalibPayload::MakeKneePayloads - run " << r.run << ": fit " << r.knee << ", sigma " << r.sigma
                << ", " << r.events << " events -> median knee " << knee_med << std::endl;
    }
  }
  std::cout << "RhoUECalibPayload::MakeKneePayloads - " << rows.size() << " runs (" << nsub << " with the median knee "
            << knee_med << ") -> " << out_dir << std::endl;
  return 0;
}

TH1 *RhoUECalibPayload::KneeHisto(const std::string &sums_file)
{
  auto f = open_file(sums_file);
  if (!f)
  {
    return nullptr;
  }
  std::map<int, double> knee;
  for (auto *obj : *f->GetListOfKeys())
  {
    const std::string name = obj->GetName();
    if (name.rfind("knee_", 0) != 0)
    {
      continue;
    }
    auto *h = f->Get<TH1>(name.c_str());
    if (h->GetBinContent(2) > 0)
    {
      knee[std::stoi(name.substr(5))] = h->GetBinContent(1) / h->GetBinContent(2);
    }
  }
  if (knee.empty())
  {
    return nullptr;
  }
  const int lo = knee.begin()->first;
  const int hi = knee.rbegin()->first + 1;
  auto *h = new TH1D("mbdq_knee", "MBD charge knee per run;run;Q_{knee}", hi - lo, lo, hi);
  h->SetDirectory(nullptr);
  for (const auto &[run, q] : knee)
  {
    h->SetBinContent(h->FindBin(run), q);
  }
  return h;
}

int RhoUECalibPayload::MakeTowerShapePayload(const std::string &sums_file, const std::string &payload,
                                             const double max_half_rms, const double min_entries)
{
  auto f = open_file(sums_file);
  if (!f)
  {
    return 1;
  }
  auto *hz = f->Get<TH1>("zvtx_bins");
  auto *hq = f->Get<TH1>("mbdq_bins");
  TH1 *hknee = KneeHisto(sums_file);
  if (!hz || !hq || !hknee)
  {
    std::cout << "RhoUECalibPayload::MakeTowerShapePayload - binning or knee missing in " << sums_file << std::endl;
    return 1;
  }
  constexpr int NE = RhoUECalibBase::NETA;
  constexpr int NP = RhoUECalibBase::NPHI;
  const int nc = hz->GetNbinsX() * hq->GetNbinsX();

  CDBHistos cdb(payload);
  auto *hfmt = new TH1D("format_version", "payload format", 1, 0, 1);
  hfmt->SetBinContent(1, 2);
  auto *hzb = new TH1D("zvtx_bins", "z_{vtx} bins;z_{vtx} [cm]", hz->GetNbinsX(), hz->GetXaxis()->GetXbins()->GetArray());
  auto *hqb = new TH1D("mbdq_bins", "Q' bins;Q' = Q/Q_{knee}", hq->GetNbinsX(), hq->GetXaxis()->GetXbins()->GetArray());
  for (TH1 *h : {static_cast<TH1 *>(hfmt), static_cast<TH1 *>(hzb), static_cast<TH1 *>(hqb), hknee})
  {
    h->SetDirectory(nullptr);
    cdb.registerHisto(h);
  }
  for (int layer = 0; layer < RhoUECalibBase::NLAYER; ++layer)
  {
    const std::string lname = RhoUECalibBase::layer_names().at(layer);
    std::array<TH3F *, 2> sum{};
    std::array<TH3F *, 2> cnt{};
    for (int half = 0; half < 2; ++half)
    {
      sum.at(half) = f->Get<TH3F>(("sum_" + lname + "_h" + std::to_string(half)).c_str());
      cnt.at(half) = f->Get<TH3F>(("cnt_" + lname + "_h" + std::to_string(half)).c_str());
      if (!sum.at(half) || !cnt.at(half))
      {
        std::cout << "RhoUECalibPayload::MakeTowerShapePayload - sums of " << lname << " missing" << std::endl;
        return 1;
      }
    }
    auto *hs = new TH3F(("s_" + lname).c_str(), ("s, " + lname + ";i#eta;i#phi;class = izbin * n_{Q'} + iQ'").c_str(),
                        NE, 0, NE, NP, 0, NP, nc, 0, nc);
    hs->SetDirectory(nullptr);
    long n_cells = 0;
    long n_unset = 0;
    long n_noisy = 0;
    for (int c = 0; c < nc; ++c)
    {
      // per tower: mean E_T of both halves (h = 2) and of each half
      std::array<std::vector<double>, 3> mean;
      std::array<double, 3> layer_mean{};
      for (int h = 0; h < 3; ++h)
      {
        mean.at(h).assign(NE * NP, -1);
        double tot = 0;
        int ng = 0;
        for (int k = 0; k < NE * NP; ++k)
        {
          const int bin = sum[0]->GetBin((k / NP) + 1, (k % NP) + 1, c + 1);
          const double s = (h == 2) ? sum[0]->GetBinContent(bin) + sum[1]->GetBinContent(bin) : sum.at(h)->GetBinContent(bin);
          const double n = (h == 2) ? cnt[0]->GetBinContent(bin) + cnt[1]->GetBinContent(bin) : cnt.at(h)->GetBinContent(bin);
          if (n >= ((h == 2) ? min_entries : 0.5 * min_entries))
          {
            mean.at(h).at(k) = s / n;
            tot += s / n;
            ++ng;
          }
        }
        layer_mean.at(h) = ng > 0 ? tot / ng : 0;
      }
      // precision of this (layer, class) map from its two halves
      double d1 = 0;
      double d2 = 0;
      int nd = 0;
      for (int k = 0; k < NE * NP; ++k)
      {
        if (mean[0].at(k) < 0 || mean[1].at(k) < 0 || mean[2].at(k) <= 0 || layer_mean[0] <= 0 || layer_mean[1] <= 0 || layer_mean[2] <= 0)
        {
          continue;
        }
        const double s = mean[2].at(k) / layer_mean[2];
        const double d = ((mean[0].at(k) / layer_mean[0]) - (mean[1].at(k) / layer_mean[1])) / s;
        d1 += d;
        d2 += d * d;
        ++nd;
      }
      const double rms = nd > 1 ? std::sqrt(std::max(0., (d2 / nd) - ((d1 / nd) * (d1 / nd)))) : 1e9;
      const bool noisy = !(layer_mean[2] > 0) || rms > max_half_rms;
      n_noisy += (noisy && nd > 1) ? 1 : 0;
      for (int k = 0; k < NE * NP; ++k)
      {
        double s = 1;
        if (!noisy && mean[2].at(k) > 0 && mean[2].at(k) / layer_mean[2] > 0.05)
        {
          s = mean[2].at(k) / layer_mean[2];
        }
        else
        {
          ++n_unset;
        }
        ++n_cells;
        hs->SetBinContent((k / NP) + 1, (k % NP) + 1, c + 1, s);
      }
    }
    cdb.registerHisto(hs);
    std::cout << "RhoUECalibPayload::MakeTowerShapePayload - " << lname << ": " << n_unset << " of " << n_cells
              << " (class, tower) cells s = 1, " << n_noisy << " classes set to s = 1 (halves differ by > " << max_half_rms << ")" << std::endl;
  }
  cdb.WriteCDBHistos();
  std::cout << "RhoUECalibPayload::MakeTowerShapePayload - " << payload << std::endl;
  return 0;
}

int RhoUECalibPayload::MakeRhoUEPayload(const std::string &sums_file, const std::string &payload, const double min_events)
{
  auto f = open_file(sums_file);
  if (!f)
  {
    return 1;
  }
  auto *hz = f->Get<TH1>("zvtx_bins");
  auto *hq = f->Get<TH1>("mbdq_bins");
  auto *zsum = f->Get<TH2D>("zsum");
  auto *qsum = f->Get<TH2D>("qsum");
  auto *nsum = f->Get<TH2D>("nsum");
  TH1 *hknee = KneeHisto(sums_file);
  if (!hz || !hq || !zsum || !qsum || !nsum || !hknee)
  {
    std::cout << "RhoUECalibPayload::MakeRhoUEPayload - binning, class sums or knee missing in " << sums_file << std::endl;
    return 1;
  }
  const int nz = hz->GetNbinsX();
  const int nq = hq->GetNbinsX();
  constexpr int NE = RhoUECalibBase::NETA;

  // interpolation nodes: event-weighted mean z per z bin (all Q'), mean Q' per Q' bin (all z)
  auto *hzb = new TH1D("zvtx_bins", "z_{vtx} bins, content = event-weighted centre;z_{vtx} [cm]", nz, hz->GetXaxis()->GetXbins()->GetArray());
  auto *hqb = new TH1D("mbdq_bins", "Q' bins, content = event-weighted centre;Q' = Q/Q_{knee}", nq, hq->GetXaxis()->GetXbins()->GetArray());
  for (int iz = 1; iz <= nz; ++iz)
  {
    double s = 0;
    double n = 0;
    for (int iq = 1; iq <= nq; ++iq)
    {
      s += zsum->GetBinContent(iz, iq);
      n += nsum->GetBinContent(iz, iq);
    }
    hzb->SetBinContent(iz, n > 0 ? s / n : hzb->GetXaxis()->GetBinCenter(iz));
  }
  for (int iq = 1; iq <= nq; ++iq)
  {
    double s = 0;
    double n = 0;
    for (int iz = 1; iz <= nz; ++iz)
    {
      s += qsum->GetBinContent(iz, iq);
      n += nsum->GetBinContent(iz, iq);
    }
    hqb->SetBinContent(iq, n > 0 ? s / n : hqb->GetXaxis()->GetBinCenter(iq));
  }

  CDBHistos cdb(payload);
  auto *hfmt = new TH1D("format_version", "payload format", 1, 0, 1);
  hfmt->SetBinContent(1, 2);
  for (TH1 *h : {static_cast<TH1 *>(hfmt), static_cast<TH1 *>(hzb), static_cast<TH1 *>(hqb), hknee})
  {
    h->SetDirectory(nullptr);
    cdb.registerHisto(h);
  }
  std::vector<double> eta_edges(NE + 1);
  for (int i = 0; i <= NE; ++i)
  {
    eta_edges[i] = i;
  }
  for (int layer = 0; layer < RhoUECalibBase::NLAYER; ++layer)
  {
    const std::string lname = RhoUECalibBase::layer_names().at(layer);
    std::array<TH3D *, RhoUECalibMaker::NSUM> S{};
    for (int k = 0; k < RhoUECalibMaker::NSUM; ++k)
    {
      S.at(k) = f->Get<TH3D>((RhoUECalibMaker::sum_names().at(k) + "_" + lname).c_str());
      if (!S.at(k))
      {
        std::cout << "RhoUECalibPayload::MakeRhoUEPayload - sums of " << lname << " missing" << std::endl;
        return 1;
      }
    }
    std::array<TH3F *, 3> out{};
    const std::array<std::string, 3> what = {"w", "alpha", "beta"};
    for (int i = 0; i < 3; ++i)
    {
      out.at(i) = new TH3F((what.at(i) + "_" + lname).c_str(), (what.at(i) + ", " + lname + ";i#eta;z_{vtx} [cm];Q'").c_str(), NE, eta_edges.data(),
                           nz, hz->GetXaxis()->GetXbins()->GetArray(), nq, hq->GetXaxis()->GetXbins()->GetArray());
      out.at(i)->SetDirectory(nullptr);
    }
    long nfit = 0;
    long nfill = 0;
    for (int ie = 0; ie < NE; ++ie)
    {
      std::vector<std::array<double, 3>> val(static_cast<size_t>(nz) * nq, {NAN, NAN, NAN});
      for (int iz = 0; iz < nz; ++iz)
      {
        for (int iq = 0; iq < nq; ++iq)
        {
          const int bin = S[0]->GetBin(ie + 1, iz + 1, iq + 1);
          const double sw = S[0]->GetBinContent(bin);
          const double sx = S[1]->GetBinContent(bin);
          const double sy = S[2]->GetBinContent(bin);
          const double sxx = S[3]->GetBinContent(bin);
          const double sxy = S[4]->GetBinContent(bin);
          const double n = S[6]->GetBinContent(bin);
          const double den = (sw * sxx) - (sx * sx);
          if (n < min_events || sw <= 0 || den <= 0)
          {
            continue;
          }
          const double alpha = ((sw * sxy) - (sx * sy)) / den;
          val.at((iz * nq) + iq) = {sx != 0 ? sy / sx : NAN, alpha, (sy - (alpha * sx)) / sw};
          ++nfit;
        }
      }
      // bins without a fit take the nearest fit of the same z row in Q', else the nearest fit of
      // the same Q' column in z; sources are the original fits only, never a filled bin
      const std::vector<std::array<double, 3>> fit = val;
      for (int iz = 0; iz < nz; ++iz)
      {
        for (int iq = 0; iq < nq; ++iq)
        {
          auto &v = val.at((iz * nq) + iq);
          if (!std::isnan(v[1]))
          {
            continue;
          }
          for (int d = 1; d < nq && std::isnan(v[1]); ++d)
          {
            for (const int j : {iq - d, iq + d})
            {
              if (j >= 0 && j < nq && !std::isnan(fit.at((iz * nq) + j)[1]))
              {
                v = fit.at((iz * nq) + j);
                break;
              }
            }
          }
          for (int d = 1; d < nz && std::isnan(v[1]); ++d)
          {
            for (const int j : {iz - d, iz + d})
            {
              if (j >= 0 && j < nz && !std::isnan(fit.at((j * nq) + iq)[1]))
              {
                v = fit.at((j * nq) + iq);
                break;
              }
            }
          }
          ++nfill;
        }
      }
      for (int iz = 0; iz < nz; ++iz)
      {
        for (int iq = 0; iq < nq; ++iq)
        {
          const auto &v = val.at((iz * nq) + iq);
          out[0]->SetBinContent(ie + 1, iz + 1, iq + 1, std::isnan(v[0]) ? 1 : v[0]);
          out[1]->SetBinContent(ie + 1, iz + 1, iq + 1, std::isnan(v[1]) ? 1 : v[1]);
          out[2]->SetBinContent(ie + 1, iz + 1, iq + 1, std::isnan(v[2]) ? 0 : v[2]);
        }
      }
    }
    for (auto *h : out)
    {
      cdb.registerHisto(h);
    }
    std::cout << "RhoUECalibPayload::MakeRhoUEPayload - " << lname << ": " << nfit << " (ieta, z, Q') bins fitted, "
              << nfill << " filled from neighbours" << std::endl;
  }
  cdb.WriteCDBHistos();
  std::cout << "RhoUECalibPayload::MakeRhoUEPayload - " << payload << std::endl;
  return 0;
}
