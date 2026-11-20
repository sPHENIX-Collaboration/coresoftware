#include "TpcDistortionQA.h"

#include <TDirectory.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3.h>
#include <TNamed.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace TpcCalibrationQADetail
{
  Fit TpcDistortionQA::fit(TH1 *h, int minEntries, bool background)
  {
    Fit f;
    f.n = h->Integral(1, h->GetNbinsX());
    if (f.n < minEntries)
    {
      return f;
    }
    double probs[] = {.16, .5, .84};
    double qs[3];
    h->GetQuantiles(3, qs, probs);
    const double width = h->GetXaxis()->GetBinWidth(h->GetMaximumBin());
    double sigma = std::max((qs[2] - qs[0]) / 2., width);
    double mean = qs[1];
    if (background)
    {
      // Locate the peak core instead of allowing broad tails to set its width.
      const int peak = h->GetMaximumBin();
      std::vector<double> contents;
      for (int b = 1; b <= h->GetNbinsX(); ++b)
      {
        contents.push_back(h->GetBinContent(b));
      }
      std::sort(contents.begin(), contents.end());
      const double pedestal = contents[contents.size() / 4];
      const double half = pedestal + .5 * (h->GetBinContent(peak) - pedestal);
      int left = peak;
      int right = peak;
      while (left > 1 && h->GetBinContent(left) > half)
      {
        --left;
      }
      while (right < h->GetNbinsX() && h->GetBinContent(right) > half)
      {
        ++right;
      }
      if (left < peak && right > peak && h->GetBinContent(left) <= half &&
          h->GetBinContent(right) <= half)
      {
        const auto crossing = [&](int a, int b)
        {
          const double ya = h->GetBinContent(a);
          const double yb = h->GetBinContent(b);
          return h->GetBinCenter(a) +
                 (half - ya) * (h->GetBinCenter(b) - h->GetBinCenter(a)) / (yb - ya);
        };
        const double low = crossing(left, left + 1);
        const double high = crossing(right - 1, right);
        sigma = std::max(width, (high - low) / 2.354820045);
        mean = h->GetBinCenter(peak);
      }
    }
    const double axisLow = h->GetXaxis()->GetXmin();
    const double axisHigh = h->GetXaxis()->GetXmax();
    if (!(sigma > 0) || !std::isfinite(mean))
    {
      f.status = 2;
      return f;
    }
    for (int pass = 0; pass < 2; ++pass)
    {
      const double half = std::max((pass ? 2. : 2.5) * sigma, (background ? 3.5 : 3.) * width);
      int lo = std::max(1, h->GetXaxis()->FindFixBin(mean - half));
      int hi = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(mean + half));
      f.low = h->GetXaxis()->GetBinLowEdge(lo);
      f.high = h->GetXaxis()->GetBinUpEdge(hi);
      if (hi - lo + 1 < 6 || h->Integral(lo, hi) < minEntries)
      {
        f.status = 2;
        return f;
      }
      TF1 model("tpcdqa_gaussian", background ? "gaus(0)+[3]" : "gaus", f.low, f.high);
      const double pedestal =
          background ? std::max(0., std::min(h->GetBinContent(lo), h->GetBinContent(hi))) : 0.;
      model.SetParameters(std::max(1., h->GetMaximum() - pedestal), mean, sigma);
      if (background)
      {
        model.SetParameter(3, pedestal);
        model.SetParLimits(3, 0., std::max(1., h->GetMaximum()));
      }
      model.SetParLimits(0, 1e-9, std::max(10., 100. * h->GetMaximum()));
      model.SetParLimits(1, f.low, f.high);
      model.SetParLimits(2, .25 * width, axisHigh - axisLow);
      // Poisson likelihood, integrate across bins, retain result, no attachment.
      TFitResultPtr result = h->Fit(&model, "QSN0LRI");
      f.rootStatus = int(result);
      if (f.rootStatus != 0 || !result.Get() || !result->IsValid())
      {
        f.status = 3;
        return f;
      }
      f.covariance = result->CovMatrixStatus();
      f.ndf = result->Ndf();
      if (f.covariance < 2)
      {
        f.status = 4;
        return f;
      }
      mean = model.GetParameter(1);
      sigma = std::abs(model.GetParameter(2));
      f.mean = mean;
      f.amplitude = model.GetParameter(0);
      f.background = background ? model.GetParameter(3) : 0.;
      f.sigma = sigma;
      f.meanError = model.GetParError(1);
      f.sigmaError = model.GetParError(2);
      f.chi2ndf = result->Ndf() > 0 ? result->Chi2() / result->Ndf() : -1.;
      if (!std::isfinite(mean) || !std::isfinite(sigma) || !std::isfinite(f.meanError) ||
          !std::isfinite(f.sigmaError) || f.meanError <= 0 || f.sigmaError <= 0 ||
          sigma <= .26 * width || mean <= f.low || mean >= f.high)
      {
        f.status = 5;
        return f;
      }
      // The seed window may contain tails. Its poor chi2 must not prevent the
      // narrower second-pass core fit; apply goodness of fit to the final pass.
      if (pass == 1 && (result->Ndf() < 3 || !std::isfinite(f.chi2ndf) || f.chi2ndf > 10. ||
                        (background && model.GetParameter(0) <= 3. * model.GetParError(0))))
      {
        f.status = 6;
        return f;
      }
    }
    f.status = 0;
    return f;
  }

  void slices(TH2 *h, const std::string &key, const std::string &title, Curves &curves,
              const Sample &sample, const Config &cfg, TDirectory *dir, std::ostream &csv, ObjectIndex &outputIndex)
  {
    Curve mean = makeCurve(title, h->GetXaxis()->GetTitle(), "Gaussian core mean DCA_{xy} [cm]");
    Curve sigma = makeCurve(title, mean.xTitle, "Gaussian core #sigma DCA_{xy} [cm]", "ratio");
    mean.edges = sigma.edges = edges(h->GetXaxis());
    const auto target = sliceLocation(key);
    outputIndex[key] = target.folder + "/" + target.name;
    auto *checks = folder(dir, "diagnostics/" + target.folder + "/" + target.name);
    auto *sliceDir = cfg.saveSlices ? folder(checks, "slices") : nullptr;
    TH1D status((key + "_fit_status").c_str(), "Fit status (0=accepted);slice coordinate;status",
                h->GetNbinsX(), mean.edges.data());
    status.SetDirectory(nullptr);
    Hists diagnostics;
    for (const std::string name : {"counts", "chi2_ndf", "root_status", "covariance", "low", "high",
                                   "underflow", "overflow", "ndf", "background"})
    {
      Hist item(new TH1D(
          (std::string(key).append("_fit_").append(name)).c_str(),
          (std::string("Slice fit ").append(name).append(";slice coordinate;").append(name)).c_str(),
          h->GetNbinsX(), mean.edges.data()));
      item->SetDirectory(nullptr);
      diagnostics.emplace(name, std::move(item));
    }
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
      Hist slice(h->ProjectionY(("slice_" + std::to_string(i)).c_str(), i, i, "e"));
      slice->SetDirectory(nullptr);
      const Fit f = TpcDistortionQA::fit(slice.get(), cfg.minEntries, cfg.fitBackground);
      status.SetBinContent(i, f.status);
      diagnostics.at("counts")->SetBinContent(i, f.n);
      diagnostics.at("ndf")->SetBinContent(i, f.ndf);
      diagnostics.at("background")->SetBinContent(i, f.background);
      diagnostics.at("chi2_ndf")->SetBinContent(i, f.chi2ndf);
      diagnostics.at("root_status")->SetBinContent(i, f.rootStatus);
      diagnostics.at("covariance")->SetBinContent(i, f.covariance);
      diagnostics.at("low")->SetBinContent(i, f.low);
      diagnostics.at("high")->SetBinContent(i, f.high);
      diagnostics.at("underflow")->SetBinContent(i, slice->GetBinContent(0));
      diagnostics.at("overflow")->SetBinContent(i, slice->GetBinContent(slice->GetNbinsX() + 1));
      const double x = h->GetXaxis()->GetBinCenter(i);
      const double ex = h->GetXaxis()->GetBinWidth(i) / 2.;
      csv << sample.label << ',' << key << ',' << i << ',' << x << ',' << ex << ',' << f.n << ','
          << f.status << ',' << f.rootStatus << ',' << f.covariance << ',' << f.mean << ','
          << f.meanError << ',' << f.sigma << ',' << f.sigmaError << ',' << f.low << ',' << f.high
          << ',' << f.chi2ndf << ',' << slice->GetBinContent(0) << ','
          << slice->GetBinContent(slice->GetNbinsX() + 1) << ',' << f.ndf << ',' << f.background
          << ',' << (cfg.fitBackground ? "gaussian_background" : "gaussian") << '\n';
      if (sliceDir)
      {
        sliceDir->cd();
        slice->Write();
      }
      if (!f.status)
      {
        TF1 fit("accepted_gaussian", cfg.fitBackground ? "gaus(0)+[3]" : "gaus", f.low, f.high);
        fit.SetParameters(f.amplitude, f.mean, f.sigma);
        if (cfg.fitBackground)
        {
          fit.SetParameter(3, f.background);
        }
        if (sliceDir)
        {
          fit.Write(("fit_" + std::to_string(i)).c_str());
        }
        point(mean, x, ex, f.mean, f.meanError);
        point(sigma, x, ex, f.sigma, f.sigmaError);
      }
    }
    checks->cd();
    status.Write("fit_status");
    for (const auto &entry : diagnostics)
    {
      entry.second->Write(entry.first.c_str());
    }
    // Keep the measured DCA map even when detailed per-bin slices are disabled.
    folder(dir, target.folder)->cd();
    h->SetTitle((plain(target.name) + " (" + plain(target.folder) + ");" + h->GetXaxis()->GetTitle() +
                 ";DCA_{xy} to origin [cm]")
                    .c_str());
    h->SetOption("COLZ");
    h->SetStats(false);
    h->Write(target.name.c_str());
    if (key.find("_phi_") != std::string::npos && key.find("phiBin") == std::string::npos)
    {
      const std::string charge = ends(key, "_pos") ? "pos" : "neg";
      Hist distribution(
          h->ProjectionY(("rdca_distribution_" + charge).c_str(), 1, h->GetNbinsX(), "e"));
      distribution->SetDirectory(nullptr);
      distribution->SetOption("HIST");
      distribution->SetTitle(("DCA_{xy} distribution, " + charge + " (" + plain(target.folder) +
                              ");DCA_{xy} to origin [cm];tracks")
                                 .c_str());
      distribution->Write();
    }
    curves[key + "_mean"] = std::move(mean);
    curves[key + "_sigma"] = std::move(sigma);
  }

  void chargeCurves(TH1 *pos, TH1 *neg, const std::string &key, const std::string &title,
                    Curves &curves, double frames)
  {
    if (!sameBins(pos, neg))
    {
      throw std::runtime_error("Charge histogram binning mismatch: " + key);
    }
    Curve ratio = makeCurve(title, pos->GetXaxis()->GetTitle(), "N_{+}/N_{-}", "ratio");
    Curve asym = makeCurve(title, ratio.xTitle, "(N_{+}-N_{-})/(N_{+}+N_{-})");
    Curve plus = makeCurve(title, ratio.xTitle, "tracks / frame / axis unit", "ratio");
    Curve minus = plus;
    Curve plusCounts = makeCurve(title, ratio.xTitle, "positive tracks", "ratio");
    Curve minusCounts = makeCurve(title, ratio.xTitle, "negative tracks", "ratio");
    ratio.edges = asym.edges = plus.edges = minus.edges = edges(pos->GetXaxis());
    plusCounts.edges = minusCounts.edges = ratio.edges;
    for (int i = 1; i <= pos->GetNbinsX(); ++i)
    {
      const double p = pos->GetBinContent(i);
      const double n = neg->GetBinContent(i);
      const double ep = pos->GetBinError(i);
      const double en = neg->GetBinError(i);
      const double x = pos->GetBinCenter(i);
      const double ex = pos->GetBinWidth(i) / 2.;
      if (p < 0 || n < 0)
      {
        throw std::runtime_error("Negative counts in " + key);
      }
      point(plusCounts, x, ex, p, ep);
      point(minusCounts, x, ex, n, en);
      if (p + n > 0)
      {
        point(asym, x, ex, (p - n) / (p + n), 2. * std::hypot(n * ep, p * en) / ((p + n) * (p + n)));
      }
      if (n > 0)
      {
        point(ratio, x, ex, p / n, std::hypot(ep / n, p * en / (n * n)));
      }
      if (frames > 0)
      {
        const double norm = frames * pos->GetBinWidth(i);
        point(plus, x, ex, p / norm, ep / norm);
        point(minus, x, ex, n / norm, en / norm);
      }
    }
    curves[key + "_charge_ratio"] = ratio;
    curves[key + "_charge_asymmetry"] = asym;
    curves[key + "_pos_counts"] = std::move(plusCounts);
    curves[key + "_neg_counts"] = std::move(minusCounts);
    if (frames > 0)
    {
      curves[key + "_pos_yield"] = plus;
      curves[key + "_neg_yield"] = minus;
    }
  }

  std::pair<int, int> ptRange(const TAxis *axis, double low, double high)
  {
    int a = -1;
    int b = -1;
    for (int i = 1; i <= axis->GetNbins() + 1; ++i)
    {
      if (std::abs(axis->GetBinLowEdge(i) - low) < 1e-7)
      {
        a = i;
      }
      if (std::abs(axis->GetBinLowEdge(i) - high) < 1e-7)
      {
        b = i - 1;
      }
    }
    if (a < 1 || b < a)
    {
      throw std::runtime_error("Requested pT limits must coincide with histogram bin edges");
    }
    return {a, b};
  }

  void rebin(TH2 *h, int group)
  {
    if (group < 1 || h->GetNbinsX() % group)
    {
      throw std::runtime_error("Slice grouping must divide the source axis bins");
    }
    if (group > 1)
    {
      h->RebinX(group);
    }
  }

  void rebin1(TH1 *h, int group)
  {
    if (group < 1 || h->GetNbinsX() % group)
    {
      throw std::runtime_error("Charge grouping must divide the source axis bins");
    }
    if (group > 1)
    {
      h->Rebin(group);
    }
  }

  void phiQoverPt(TH3 *pos, TH3 *neg, const std::string &side, const Sample &sample,
                  const Config &cfg, TDirectory *dir, Curves &curves, std::ostream &csv, ObjectIndex &outputIndex)
  {
    if (!sameBins(pos, neg) || pos->GetNbinsX() % cfg.phiGroup)
    {
      throw std::runtime_error("Incompatible phi-q/pT source bins");
    }
    std::vector<double> inverse;
    for (int j = 1; j <= pos->GetNbinsY() + 1; ++j)
    {
      const double pt = pos->GetYaxis()->GetBinLowEdge(j);
      if (pt >= .2 - 1e-7)
      {
        inverse.push_back(1. / pt);
      }
    }
    std::vector<double> qedges;
    qedges.reserve(2 * inverse.size());
    for (double v : inverse)
    {
      qedges.push_back(-v);
    }
    qedges.insert(qedges.end(), inverse.rbegin(), inverse.rend());
    const auto dcaEdges = edges(pos->GetZaxis());
    for (int lo = 1; lo <= pos->GetNbinsX(); lo += cfg.phiGroup)
    {
      const int hi = lo + cfg.phiGroup - 1;
      const std::string key =
          "new_" + side + "_phiBin" + std::to_string((lo - 1) / cfg.phiGroup + 1) + "_qoverpt";
      TH2D counts((key + "_counts").c_str(), ";q/p_{T} [c/GeV];DCA_{xy} [cm]", int(qedges.size()) - 1,
                  qedges.data(), int(dcaEdges.size()) - 1, dcaEdges.data());
      counts.SetDirectory(nullptr);
      for (const int sign : {-1, 1})
      {
        auto *source = sign > 0 ? pos : neg;
        for (int j = 1; j <= source->GetNbinsY(); ++j)
        {
          const double pl = source->GetYaxis()->GetBinLowEdge(j);
          const double ph = source->GetYaxis()->GetBinUpEdge(j);
          if (pl < .2 - 1e-7)
          {
            continue;
          }
          const int qbin = counts.GetXaxis()->FindFixBin(sign * .5 * (1. / pl + 1. / ph));
          Hist dca(source->ProjectionZ("phi_qpt_slice", lo, hi, j, j, "e"));
          dca->SetDirectory(nullptr);
          for (int k = 0; k <= dca->GetNbinsX() + 1; ++k)
          {
            counts.SetBinContent(qbin, k, dca->GetBinContent(k));
            counts.SetBinError(qbin, k, dca->GetBinError(k));
          }
        }
      }
      std::ostringstream title;
      title << sample.label << ": TPC-only hardware " << side << ", "
            << pos->GetXaxis()->GetBinLowEdge(lo) << "<=phi<" << pos->GetXaxis()->GetBinUpEdge(hi)
            << ", all stored pT>0.2";
      slices(&counts, key, title.str(), curves, sample, cfg, dir, csv, outputIndex);
    }
  }

  Curves TpcDistortionQA::analyze(Hists &hs, const Sample &sample, const Config &cfg, TDirectory *dir,
                                  std::ostream &csv, std::ostream &trends, ObjectIndex &outputIndex)
  {
    Curves curves;
    const auto get = [&](const std::string &name) -> TH1 *
    {
      auto i = hs.find(name);
      return i == hs.end() ? nullptr : i->second.get();
    };
    const double frames = get("nrecotracks1d") ? get("nrecotracks1d")->GetEntries() : 0.;
    folder(dir, "metadata")->cd();
    TNamed exposure("exposure",
                    ("QA frames=" + std::to_string(frames) + "; input-list=" + sample.list).c_str());
    exposure.Write();
    for (const std::string n : {"tpc_only_crossing", "tpc_only_status"})
    {
      if (get(n))
      {
        folder(dir, "diagnostics/selection")->cd();
        get(n)->Write(n == "tpc_only_crossing" ? "assigned_crossing" : "track_selection_status");
      }
    }
    // Legacy profiles always remain separate from hardware-side, selected new QA.
    for (const std::string side : {"south", "north"})
    {
      const std::string legacyTitle = sample.label + ": legacy PCA-z " +
                                      (side == "south" ? "<=0" : ">0") + ", no new pT/cluster cut";
      for (const std::string charge : {"pos", "neg"})
      {
        auto *h = dynamic_cast<TH2 *>(
            get(std::string("dcaxyorigin_phi_").append(side).append("_").append(charge)));
        if (!h)
        {
          continue;
        }
        Hist copy(static_cast<TH1 *>(h->Clone()));
        copy->SetDirectory(nullptr);
        auto *h2 = static_cast<TH2 *>(copy.get());
        rebin(h2, cfg.phiGroup);
        slices(h2, std::string("legacy_pcaz_").append(side).append("_phi_").append(charge),
               legacyTitle, curves, sample, cfg, dir, csv, outputIndex);
      }
      const std::string newTitle =
          sample.label + ": TPC-only hardware " + side + ", N_{TPC}>=31, all assigned crossings";
      for (const std::string axis : {"phi", "eta"})
      {
        for (const std::string charge : {"pos", "neg"})
        {
          auto *h = dynamic_cast<TH3 *>(get(std::string("tpc_only_")
                                                .append(side)
                                                .append("_")
                                                .append(charge)
                                                .append("_rdca_")
                                                .append(axis)
                                                .append("_pt")));
          if (!h)
          {
            continue;
          }
          const auto bins = ptRange(h->GetYaxis(), cfg.ptMin, cfg.ptMax);
          h->GetYaxis()->SetRange(bins.first, bins.second);
          std::unique_ptr<TH2> projected(static_cast<TH2 *>(h->Project3D("zxe")));
          projected->SetDirectory(nullptr);
          h->GetYaxis()->SetRange(0, 0);
          rebin(projected.get(), axis == "phi" ? cfg.phiGroup : cfg.etaGroup);
          slices(projected.get(),
                 std::string("new_").append(side).append("_").append(axis).append("_").append(charge),
                 newTitle + ", " + std::to_string(cfg.ptMin) + "<=pT<" + std::to_string(cfg.ptMax),
                 curves, sample, cfg, dir, csv, outputIndex);
        }
      }
      if (auto *h = dynamic_cast<TH2 *>(get("tpc_only_" + side + "_rdca_qoverpt")))
      {
        // No extra pT cut is possible here: the source already integrates pT>0.2.
        slices(h, "new_" + side + "_qoverpt", newTitle + ", pT>0.2 (not the selected pT window)",
               curves, sample, cfg, dir, csv, outputIndex);
      }
      if (cfg.phiQpt)
      {
        auto *p = dynamic_cast<TH3 *>(get("tpc_only_" + side + "_pos_rdca_phi_pt"));
        auto *n = dynamic_cast<TH3 *>(get("tpc_only_" + side + "_neg_rdca_phi_pt"));
        if (p && n)
        {
          phiQoverPt(p, n, side, sample, cfg, dir, curves, csv, outputIndex);
        }
      }
      auto *p3 = dynamic_cast<TH3 *>(get("tpc_only_" + side + "_pos_pt_phi_eta"));
      auto *n3 = dynamic_cast<TH3 *>(get("tpc_only_" + side + "_neg_pt_phi_eta"));
      if (p3 && n3)
      {
        if (!sameBins(p3, n3))
        {
          throw std::runtime_error("New charge source bins differ");
        }
        const auto bins = ptRange(p3->GetZaxis(), cfg.ptMin, cfg.ptMax);
        for (const std::string axis : {"phi", "eta"})
        {
          Hist p(axis == "phi"
                     ? p3->ProjectionX("charge_p", 1, p3->GetNbinsY(), bins.first, bins.second, "e")
                     : p3->ProjectionY("charge_p", 1, p3->GetNbinsX(), bins.first, bins.second, "e"));
          p->SetDirectory(nullptr);
          Hist n(axis == "phi"
                     ? n3->ProjectionX("charge_n", 1, n3->GetNbinsY(), bins.first, bins.second, "e")
                     : n3->ProjectionY("charge_n", 1, n3->GetNbinsX(), bins.first, bins.second, "e"));
          n->SetDirectory(nullptr);
          rebin1(p.get(), axis == "phi" ? cfg.phiGroup : cfg.etaGroup);
          rebin1(n.get(), axis == "phi" ? cfg.phiGroup : cfg.etaGroup);
          chargeCurves(
              p.get(), n.get(), std::string("new_").append(side).append("_").append(axis),
              newTitle + ", selected pT window, charge projections within |eta|<1.2, |phi|<pi",
              curves, frames);
        }
      }
    }
    auto *lp = dynamic_cast<TH2 *>(get("nrecotracks_pos"));
    auto *ln = dynamic_cast<TH2 *>(get("nrecotracks_neg"));
    if (lp && ln)
    {
      if (!sameBins(lp, ln))
      {
        throw std::runtime_error("Legacy charge source bins differ");
      }
      for (const std::string axis : {"phi", "eta"})
      {
        Hist p(axis == "eta" ? lp->ProjectionX("legacy_p", 1, lp->GetNbinsY(), "e")
                             : lp->ProjectionY("legacy_p", 1, lp->GetNbinsX(), "e"));
        p->SetDirectory(nullptr);
        Hist n(axis == "eta" ? ln->ProjectionX("legacy_n", 1, ln->GetNbinsY(), "e")
                             : ln->ProjectionY("legacy_n", 1, ln->GetNbinsX(), "e"));
        n->SetDirectory(nullptr);
        rebin1(p.get(), axis == "phi" ? cfg.phiGroup : cfg.etaGroup);
        rebin1(n.get(), axis == "phi" ? cfg.phiGroup : cfg.etaGroup);
        chargeCurves(p.get(), n.get(), "legacy_" + axis,
                     sample.label + ": legacy track map, pT>1, within stored eta/phi axes", curves,
                     frames);
      }
    }
    if (get("pt_pos") && get("pt_neg"))
    {
      chargeCurves(get("pt_pos"), get("pt_neg"), "legacy_pt",
                   sample.label + ": legacy track map, all pT", curves, frames);
    }
    // Charge-even/odd are HALF sum/difference of separately fitted core means.
    const Curves initial = curves;
    for (const auto &entry : initial)
    {
      const auto &key = entry.first;
      if (key.size() < 9 || key.substr(key.size() - 9) != "_pos_mean")
      {
        continue;
      }
      const auto base = key.substr(0, key.size() - 9);
      auto neg = initial.find(base + "_neg_mean");
      if (neg == initial.end())
      {
        continue;
      }
      for (const std::string op : {"even", "odd"})
      {
        Curve c = combine(entry.second, neg->second, op);
        c.yTitle = op == "even" ? "(#mu_{+}+#mu_{-})/2 [cm]" : "(#mu_{+}-#mu_{-})/2 [cm]";
        curves[std::string(base).append("_charge_").append(op)] = std::move(c);
      }
    }
    dir->cd();
    for (const auto &entry : curves)
    {
      auto g = graph(entry.second, entry.first);
      writeCurve(entry.second, entry.first, dir, outputIndex);
      // Harmonics are diagnostics, not a conversion to a space-charge constant.
      if (entry.first.find("phi_charge_") != std::string::npos &&
          (entry.first.find("_even") != std::string::npos ||
           entry.first.find("_odd") != std::string::npos) &&
          g->GetN() >= 8)
      {
        TF1 harmonics("harmonics", "[0]+[1]*sin(x)+[2]*cos(x)+[3]*sin(2*x)+[4]*cos(2*x)",
                      -std::acos(-1.), std::acos(-1.));
        const auto fit = g->Fit(&harmonics, "QSN0EX0R");
        if (int(fit) == 0 && fit.Get() && fit->CovMatrixStatus() >= 2)
        {
          const auto address = location(entry.first);
          folder(dir, "diagnostics/" + address.folder + "/calibration_fits")->cd();
          fit->Write((address.name + "_harmonic_result").c_str());
          harmonics.Write((address.name + "_harmonics").c_str());
          for (int p = 0; p < 5; ++p)
          {
            trends << sample.label << ',' << entry.first << ",harmonic_" << p << ','
                   << harmonics.GetParameter(p) << ',' << harmonics.GetParError(p) << '\n';
          }
        }
      }
      if (entry.first.find("qoverpt_mean") != std::string::npos && g->GetN() >= 6)
      {
        // Fit each charge separately; don't bridge the unpopulated high-pT gap.
        for (const int sign : {-1, 1})
        {
          TGraphErrors subset;
          for (size_t i = 0; i < entry.second.x.size(); ++i)
          {
            if (sign * entry.second.x[i] > 0)
            {
              const int j = subset.GetN();
              subset.SetPoint(j, entry.second.x[i], entry.second.y[i]);
              subset.SetPointError(j, 0, entry.second.ey[i]);
            }
          }
          if (subset.GetN() < 4)
          {
            continue;
          }
          TF1 line("line", "pol1", sign < 0 ? -5 : 0, sign < 0 ? 0 : 5);
          const auto fit = subset.Fit(&line, "QSN0R");
          if (int(fit) == 0 && fit.Get() && fit->CovMatrixStatus() >= 2)
          {
            const std::string name = entry.first + (sign < 0 ? "_neg" : "_pos");
            const auto address = location(entry.first);
            folder(dir, "diagnostics/" + address.folder + "/calibration_fits")->cd();
            const std::string publicName = address.name + (sign < 0 ? "_neg" : "_pos");
            outputIndex[name] =
                "diagnostics/" + address.folder + "/calibration_fits/" + publicName + "_line_result";
            fit->Write((publicName + "_line_result").c_str());
            line.Write((publicName + "_line").c_str());
            for (int p = 0; p < 2; ++p)
            {
              trends << sample.label << ',' << name << ',' << (p ? "slope" : "intercept") << ','
                     << line.GetParameter(p) << ',' << line.GetParError(p) << '\n';
            }
          }
        }
      }
    }
    return curves;
  }
}  // namespace TpcCalibrationQADetail
