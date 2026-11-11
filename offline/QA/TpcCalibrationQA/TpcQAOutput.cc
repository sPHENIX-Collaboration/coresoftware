#include "TpcQAInternal.h"

namespace TpcCalibrationQADetail
{
  std::map<std::string, std::string> outputIndex;

  Curve makeCurve(std::string title, std::string xTitle, std::string yTitle,
                  std::string comparison)
  {
    return {std::move(title), std::move(xTitle), std::move(yTitle), std::move(comparison), {}, {}, {}, {}, {}};
  }
  void point(Curve &c, double x, double ex, double y, double ey)
  {
    c.x.push_back(x);
    c.ex.push_back(ex);
    c.y.push_back(y);
    c.ey.push_back(ey);
  }

  std::unique_ptr<TGraphErrors> graph(const Curve &c, const std::string &name)
  {
    auto g = std::make_unique<TGraphErrors>(int(c.x.size()), c.x.data(), c.y.data(), c.ex.data(),
                                            c.ey.data());
    g->SetName(name.c_str());
    g->SetTitle((c.title + ";" + c.xTitle + ";" + c.yTitle).c_str());
    return g;
  }

  std::vector<double> edges(const TAxis *axis)
  {
    std::vector<double> e;
    for (int i = 1; i <= axis->GetNbins() + 1; ++i)
    {
      e.push_back(axis->GetBinLowEdge(i));
    }
    return e;
  }

  TDirectory *folder(TDirectory *root, const std::string &path)
  {
    std::istringstream parts(path);
    std::string part;
    auto *dir = root;
    while (std::getline(parts, part, '/'))
    {
      if (part.empty())
      {
        continue;
      }
      auto *next = dir->GetDirectory(part.c_str());
      if (!next)
      {
        next = dir->mkdir(part.c_str());
      }
      if (!next)
      {
        throw std::runtime_error("Cannot create output folder: " + path);
      }
      dir = next;
    }
    return dir;
  }

  bool ends(const std::string &s, const std::string &suffix)
  {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
  }

  std::string plain(std::string s)
  {
    std::replace(s.begin(), s.end(), '_', ' ');
    return s;
  }

  Location location(std::string key, const std::string &comparison)
  {
    std::string suffix;
    if (ends(key, "_ratio_to_reference"))
    {
      key.resize(key.size() - 19);
      suffix = "_ratio_to_reference";
    }
    else if (ends(key, "_vs_reference"))
    {
      key.resize(key.size() - 13);
      suffix = comparison == "ratio" ? "_ratio_to_reference" : "_difference_from_reference";
    }
    if (key.rfind("drift_", 0) == 0)
    {
      const std::map<std::string, std::string> names = {
          {"drift_t0_bins", "t0_vs_side"},
          {"drift_cm_bins", "cm_time_vs_side"},
          {"drift_data_velocity", "drift_velocity_data_vs_side"},
          {"drift_phg_reference_velocity", "drift_velocity_phgarfield_reference_vs_side"},
          {"drift_phg_configured_t0", "t0_phgarfield_configured_vs_side"},
          {"drift_data_over_reference", "velocity_ratio_data_over_phgarfield_reference_vs_side"},
          {"drift_phg_edge_velocity", "drift_velocity_phgarfield_edges_vs_side"},
          {"drift_phg_baseline_velocity", "drift_velocity_phgarfield_baseline_vs_side"},
          {"drift_phg_adc_weighted_velocity", "drift_velocity_phgarfield_adc_weighted_vs_side"},
          {"drift_phg_over_data", "velocity_ratio_phgarfield_edges_over_data_vs_side"},
          {"drift_data_over_baseline", "velocity_ratio_data_over_phgarfield_baseline_vs_side"},
          {"drift_phg_cm_z", "phgarfield_cm_edge_z_vs_side"},
          {"drift_phg_pad_z", "phgarfield_pad_edge_z_vs_side"}};
      auto found = names.find(key);
      std::string name = found == names.end() ? key.substr(6) : found->second;
      for (const int side : {0, 1})
      {
        const std::string old = "_profile_side" + std::to_string(side);
        if (ends(name, old))
        {
          name = name.substr(0, name.size() - old.size()) + "_edge_profile_" +
                 (side == 0 ? "south" : "north");
        }
      }
      return {suffix.empty() ? "drift" : "drift/reference_comparison", name + suffix};
    }
    std::vector<std::string> tokens;
    std::istringstream in(key);
    std::string token;
    while (std::getline(in, token, '_'))
    {
      tokens.push_back(token);
    }
    std::string selection;
    size_t start = 0;
    bool dca = false;
    if (tokens.size() >= 4 && tokens[0] == "legacy" && tokens[1] == "pcaz")
    {
      selection = tokens[2] == "south" ? "legacy_pca_z_nonpositive" : "legacy_pca_z_positive";
      start = 3;
      dca = true;
    }
    else if (tokens.size() >= 3 && tokens[0] == "new")
    {
      selection = "tpc_only_" + tokens[1];
      start = 2;
    }
    else if (tokens.size() >= 3 && tokens[0] == "legacy")
    {
      selection = tokens[1] == "pt" ? "legacy_all_tracks" : "legacy_pt_above_1";
      start = 1;
    }
    else
    {
      return {"diagnostics/unclassified", key + suffix};
    }
    if (tokens[start].rfind("phiBin", 0) == 0)
    {
      std::ostringstream number;
      number << std::setw(2) << std::setfill('0') << std::stoi(tokens[start].substr(6));
      selection += "/phi_bin_" + number.str();
      ++start;
      dca = true;
    }
    const std::string axis = tokens[start];
    const std::string last = tokens.back();
    dca = dca || last == "mean" || last == "sigma" || last == "even" || last == "odd";
    std::string name;
    if (dca)
    {
      if (last == "even" || last == "odd")
      {
        name = "rdca_charge_" + last + "_vs_" + axis;
      }
      else
      {
        name = "rdca_" + last + "_vs_" + axis;
        if (tokens.size() > start + 2 && (tokens[start + 1] == "pos" || tokens[start + 1] == "neg"))
        {
          name += "_" + tokens[start + 1];
        }
      }
    }
    else if (last == "ratio" || last == "asymmetry")
    {
      name = "charge_" + last + "_vs_" + axis;
    }
    else
    {
      name = (last == "counts" ? "track_counts_vs_" : "track_yield_per_frame_vs_") + axis + "_" +
             tokens[start + 1];
    }
    std::string path = (dca ? "rdca/" : "charge/") + selection;
    if (!suffix.empty())
    {
      path += "/reference_comparison";
    }
    // Extra representations remain available, without crowding the physics view.
    if (last == "asymmetry" || last == "yield")
    {
      path = "diagnostics/" + path;
    }
    return {path, name + suffix};
  }

  Location sliceLocation(const std::string &key)
  {
    auto out = location(key + "_mean");
    out.name.replace(0, std::string("rdca_mean_").size(), "rdca_");
    return out;
  }

  std::string writeCurve(const Curve &c, const std::string &name, TDirectory *dir, bool flat)
  {
    auto out = location(name, c.comparison);
    if (flat)
    {
      out.folder = "";
      out.name = name;
      const auto p = name.rfind("_side");
      if (p != std::string::npos)
      {
        auto metric = location(name.substr(0, p)).name;
        const auto q = metric.rfind("_vs_side");
        if (q != std::string::npos)
        {
          metric.replace(q, 8, "_vs_run");
        }
        out.name = metric + (name.substr(p + 5) == "0" ? "_south" : "_north");
      }
    }
    folder(dir, "metadata/comparisons")->cd();
    TNamed mode(name.c_str(), c.comparison.c_str());
    mode.Write();
    auto *main = folder(dir, out.folder);
    const std::string diagnostic =
        out.folder.rfind("diagnostics/", 0) == 0 ? out.folder : "diagnostics/" + out.folder;
    auto *checks = folder(dir, diagnostic);
    outputIndex[name] = out.folder + "/" + out.name;
    if (c.edges.size() < 2)
    {
      main->cd();
      auto g = graph(c, out.name);
      g->SetTitle((plain(out.name) + ";" + c.xTitle + ";" + c.yTitle).c_str());
      g->Write();
      return out.folder + "/" + out.name;
    }
    TH1D h(out.name.c_str(),
           (plain(out.name) + " (" + plain(out.folder) + ");" + c.xTitle + ";" + c.yTitle).c_str(),
           int(c.edges.size()) - 1, c.edges.data());
    TH1D mask((out.name + "_valid").c_str(), "Validity mask: 1=usable;axis;valid",
              int(c.edges.size()) - 1, c.edges.data());
    h.SetDirectory(nullptr);
    h.SetStats(false);
    h.SetOption("E1");
    mask.SetDirectory(nullptr);
    for (size_t i = 0; i < c.x.size(); ++i)
    {
      const int b = h.FindFixBin(c.x[i]);
      if (b < 1 || b > h.GetNbinsX())
      {
        continue;
      }
      h.SetBinContent(b, c.y[i]);
      h.SetBinError(b, c.ey[i]);
      mask.SetBinContent(b, 1);
    }
    main->cd();
    h.Write();
    checks->cd();
    mask.Write();
    graph(c, out.name + "_valid_points")->Write();
    folder(dir, "metadata/definitions/" + out.folder)->cd();
    TNamed definition(out.name.c_str(), c.title.c_str());
    definition.Write();
    return out.folder.empty() ? out.name : out.folder + "/" + out.name;
  }

  bool samePoint(const Curve &a, size_t i, const Curve &b, size_t j)
  {
    return std::abs(a.x[i] - b.x[j]) < 1e-9 && std::abs(a.ex[i] - b.ex[j]) < 1e-9;
  }

  Curve combine(const Curve &a, const Curve &b, const std::string &op, double ratioFloor)
  {
    Curve c = a;
    c.x.clear();
    c.ex.clear();
    c.y.clear();
    c.ey.clear();
    for (size_t i = 0; i < a.x.size(); ++i)
    {
      for (size_t j = 0; j < b.x.size(); ++j)
      {
        if (!samePoint(a, i, b, j))
        {
          continue;
        }
        double y = 0;
        double err = 0;
        if (op == "ratio" || op == "signed_ratio")
        {
          if ((op == "ratio" && b.y[j] <= 0) ||
              (op == "signed_ratio" && std::abs(b.y[j]) <= std::max(ratioFloor, 3. * b.ey[j])))
          {
            break;
          }
          y = a.y[i] / b.y[j];
          err = std::hypot(a.ey[i] / b.y[j], a.y[i] * b.ey[j] / (b.y[j] * b.y[j]));
        }
        else
        {
          const double sign = op == "even" ? 1. : -1.;
          const double scale = op == "even" || op == "odd" ? .5 : 1.;
          y = scale * (a.y[i] + sign * b.y[j]);
          err = scale * std::hypot(a.ey[i], b.ey[j]);
        }
        if (std::isfinite(y) && std::isfinite(err))
        {
          point(c, a.x[i], a.ex[i], y, err);
        }
        break;
      }
    }
    return c;
  }

  void plot(const Curves &curves, const Curves &reference, const Sample &sample,
            const std::string &refLabel, const std::string &pdf, TDirectory *dir, const Config &cfg)
  {
    TCanvas canvas("tpcdqa_canvas", "TPC distortion QA", 1100, 1100);
    const auto readable = [](std::string s)
    {
      std::replace(s.begin(), s.end(), '_', ' ');
      return s;
    };
    std::map<std::string, std::string> pdfs;
    for (const auto &entry : curves)
    {
      if (entry.second.x.empty())
      {
        continue;
      }
      const auto address = location(entry.first, entry.second.comparison);
      const auto category = address.folder.substr(0, address.folder.find('/'));
      if (category == "diagnostics")
      {
        auto referenceCurve = reference.find(entry.first);
        if (sample.label != refLabel && referenceCurve != reference.end())
        {
          auto change = combine(entry.second, referenceCurve->second, entry.second.comparison);
          writeCurve(change, entry.first + "_vs_reference", dir);
        }
        continue;
      }
      if (!hasKey(pdfs, category))
      {
        pdfs[category] = pdf.substr(0, pdf.size() - 4) + "_" + category + ".pdf";
        canvas.Print((pdfs[category] + "[").c_str());
      }
      canvas.Clear();
      const bool signedDCA = entry.first.find("_mean") != std::string::npos ||
                             entry.first.find("_charge_even") != std::string::npos ||
                             entry.first.find("_charge_odd") != std::string::npos;
      auto ref = reference.find(entry.first);
      const bool hasRef =
          sample.label != refLabel && ref != reference.end() && !ref->second.x.empty();
      int panels = 1;
      if (hasRef)
      {
        panels = signedDCA ? 3 : 2;
      }
      canvas.Divide(1, panels);
      for (int p = 1; p <= panels; ++p)
      {
        canvas.cd(p);
        gPad->SetLeftMargin(.16);
        gPad->SetRightMargin(.04);
        gPad->SetTopMargin(.12);
        gPad->SetBottomMargin(.14);
      }
      auto current = graph(entry.second, entry.first);
      std::string caption = readable(sample.label) + ": " + readable(address.name) + " (" +
                            readable(address.folder) + ")";
      if (entry.first.rfind("new_", 0) == 0 && entry.first.find("qoverpt") == std::string::npos)
      {
        std::ostringstream window;
        window << ", " << cfg.ptMin << "<=pT<" << cfg.ptMax;
        caption += window.str();
      }
      current->SetTitle((caption + ";" + entry.second.xTitle + ";" + entry.second.yTitle).c_str());
      current->SetLineColor(kBlue + 1);
      current->SetMarkerColor(kBlue + 1);
      current->SetMarkerStyle(20);
      std::unique_ptr<TGraphErrors> rg;
      std::unique_ptr<TGraphErrors> dg;
      std::unique_ptr<TGraphErrors> ratioGraph;
      canvas.cd(1);
      gPad->SetGrid();
      current->Draw("AP");
      TLegend legend(.70, .76, .97, .94);
      legend.SetTextSize(.035);
      legend.AddEntry(current.get(), sample.label.c_str(), "lp");
      if (hasRef)
      {
        rg = graph(ref->second, "reference");
        rg->SetLineColor(kRed + 1);
        rg->SetMarkerColor(kRed + 1);
        rg->SetMarkerStyle(24);
        rg->Draw("P SAME");
        legend.AddEntry(rg.get(), refLabel.c_str(), "lp");
        double low = std::numeric_limits<double>::infinity();
        double high = -low;
        for (const auto *c : {&entry.second, &ref->second})
        {
          for (size_t i = 0; i < c->y.size(); ++i)
          {
            low = std::min(low, c->y[i] - c->ey[i]);
            high = std::max(high, c->y[i] + c->ey[i]);
          }
        }
        const double margin = std::max(1e-5, .15 * (high - low));
        current->SetMinimum(low - margin);
        current->SetMaximum(high + margin);
      }
      legend.Draw();
      if (hasRef)
      {
        canvas.cd(2);
        gPad->SetGrid();
      }
      std::unique_ptr<TLine> line;
      std::unique_ptr<TLine> ratioLine;
      if (hasRef)
      {
        Curve diff = combine(entry.second, ref->second, entry.second.comparison);
        diff.title =
            readable(address.name) + ": " + readable(sample.label) + " vs " + readable(refLabel);
        diff.yTitle = entry.second.comparison == "ratio"
                          ? "sample / reference"
                          : entry.second.yTitle + " (sample - reference)";
        if (!diff.x.empty())
        {
          dg = graph(diff, entry.first + "_vs_reference");
          dg->SetMarkerStyle(20);
          if (entry.first.find("_profile_side") != std::string::npos)
          {
            double maximum = 0.;
            for (double value : diff.y)
            {
              maximum = std::max(maximum, std::abs(value));
            }
            if (maximum < .01)
            {
              dg->SetMinimum(-.01);
              dg->SetMaximum(.01);
            }
          }
          dg->Draw("AP");
          const double baseline = entry.second.comparison == "ratio" ? 1. : 0.;
          line = std::make_unique<TLine>(diff.x.front() - diff.ex.front(), baseline,
                                         diff.x.back() + diff.ex.back(), baseline);
          line->SetLineStyle(2);
          line->Draw();
          writeCurve(diff, entry.first + "_vs_reference", dir);
        }
        // Retain an empty graph/full zero mask too, when all denominators fail.
        if (diff.x.empty())
        {
          writeCurve(diff, entry.first + "_vs_reference", dir);
        }
        if (signedDCA)
        {
          Curve ratio = combine(entry.second, ref->second, "signed_ratio", cfg.ratioFloor);
          ratio.title =
              readable(address.name) + ": signed ratio, reference |mean| > max(floor, 3*error)";
          ratio.yTitle = "sample / reference (masked near zero)";
          writeCurve(ratio, entry.first + "_ratio_to_reference", dir);
          canvas.cd(3);
          gPad->SetGrid();
          if (!ratio.x.empty())
          {
            ratioGraph = graph(ratio, "ratio_display");
            ratioGraph->SetMarkerStyle(20);
            ratioGraph->Draw("AP");
            ratioLine = std::make_unique<TLine>(ratio.x.front() - ratio.ex.front(), 1.,
                                                ratio.x.back() + ratio.ex.back(), 1.);
            ratioLine->SetLineStyle(2);
            ratioLine->Draw();
          }
        }
      }
      canvas.Print(pdfs.at(category).c_str());
    }
    for (const auto &entry : pdfs)
    {
      canvas.Print((entry.second + "]").c_str());
    }
  }

  void plotDcaDistributions(TDirectory *dir, const std::string &pdf)
  {
    auto *rdca = dir->GetDirectory("rdca");
    if (!rdca)
    {
      return;
    }
    TCanvas canvas("rdca_distributions", "Measured DCA distributions", 1100, 800);
    bool opened = false;
    std::function<void(TDirectory *)> visit = [&](TDirectory *current)
    {
      TIter next(current->GetListOfKeys());
      while (auto *key = dynamic_cast<TKey *>(next()))
      {
        const std::string name = key->GetName();
        if (name == "reference_comparison")
        {
          continue;
        }
        if (auto *child = current->GetDirectory(name.c_str()))
        {
          visit(child);
          continue;
        }
        if (name.rfind("rdca_vs_", 0) != 0 && name.rfind("rdca_distribution_", 0) != 0)
        {
          continue;
        }
        std::unique_ptr<TObject> object(key->ReadObj());
        auto *hist = dynamic_cast<TH1 *>(object.get());
        if (!hist)
        {
          continue;
        }
        hist->SetDirectory(nullptr);
        if (!opened)
        {
          canvas.Print((pdf + "[").c_str());
          opened = true;
        }
        canvas.Clear();
        canvas.cd();
        canvas.SetLeftMargin(.14);
        canvas.SetRightMargin(hist->GetDimension() == 2 ? .16 : .05);
        hist->Draw(hist->GetDimension() == 2 ? "COLZ" : "HIST");
        canvas.Print(pdf.c_str());
      }
    };
    visit(rdca);
    if (opened)
    {
      canvas.Print((pdf + "]").c_str());
    }
  }
}  // namespace TpcCalibrationQADetail
