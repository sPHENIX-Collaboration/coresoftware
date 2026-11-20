#include "TpcCalibrationQA.h"
#include "TpcDistortionQA.h"
#include "TpcDriftQA.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using TpcCalibrationQADetail::combine;
using TpcCalibrationQADetail::Config;
using TpcCalibrationQADetail::Curves;
using TpcCalibrationQADetail::DriftInput;
using TpcCalibrationQADetail::ends;
using TpcCalibrationQADetail::folder;
using TpcCalibrationQADetail::graph;
using TpcCalibrationQADetail::Hists;
using TpcCalibrationQADetail::load;
using TpcCalibrationQADetail::ObjectIndex;
using TpcCalibrationQADetail::plot;
using TpcCalibrationQADetail::plotDcaDistributions;
using TpcCalibrationQADetail::point;
using TpcCalibrationQADetail::Sample;
using TpcCalibrationQADetail::TpcDistortionQA;
using TpcCalibrationQADetail::TpcDriftQA;
using TpcCalibrationQADetail::writeCurve;

int TpcCalibrationQA::run(const std::string &manifestName, const std::string &outputName,
                          const std::string &referenceName) const
{
  const char *manifest = manifestName.c_str();
  const char *output = outputName.c_str();
  const char *reference = referenceName.c_str();
  const double ptMin = m_ptMin;
  const double ptMax = m_ptMax;
  const int minEntries = m_minEntries;
  const int phiGroup = m_phiGroup;
  const int etaGroup = m_etaGroup;
  const bool saveSlices = m_saveSlices;
  const double ratioFloor = m_ratioFloor;
  const bool phiQpt = m_phiQpt;
  const double rawTimeShift = m_rawTimeShift;
  const bool makePlots = m_makePlots;
  const bool fitBackground = m_fitBackground;
  const char *prefix = m_prefix.c_str();
  const char *conditions = m_conditions.c_str();
  const char *phgReferences = m_phgReferences.c_str();
  if (!(ptMin >= .2 && ptMax > ptMin) || minEntries < 20 || phiGroup < 1 || etaGroup < 1 ||
      ratioFloor < 0 || !std::isfinite(ratioFloor) || !std::isfinite(rawTimeShift))
  {
    throw std::runtime_error("Invalid analysis configuration");
  }
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetPaperSize(20., 20.);
  TH1::AddDirectory(false);
  std::ifstream in(manifest);
  if (!in)
  {
    throw std::runtime_error("Cannot read manifest");
  }
  std::vector<Sample> samples;
  std::set<std::string> labels;
  std::string row;
  while (std::getline(in, row))
  {
    if (row.empty() || row[0] == '#')
    {
      continue;
    }
    const auto tab = row.find('\t');
    if (tab == std::string::npos)
    {
      throw std::runtime_error("Manifest needs label<TAB>filelist");
    }
    const auto second = row.find('\t', tab + 1);
    const auto third = second == std::string::npos ? std::string::npos : row.find('\t', second + 1);
    Sample s{
        row.substr(0, tab),
        row.substr(tab + 1, second == std::string::npos ? std::string::npos : second - tab - 1),
        second == std::string::npos
            ? ""
            : row.substr(second + 1,
                         third == std::string::npos ? std::string::npos : third - second - 1),
        third == std::string::npos ? std::numeric_limits<double>::quiet_NaN()
                                   : std::stod(row.substr(third + 1))};
    if (s.label.empty() || s.label == "." || s.label == ".." || s.label == "configuration" ||
        s.label == "run_conditions_json" || s.label == "drift_run_trends" ||
        s.label.find_first_not_of(
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-") !=
            std::string::npos ||
        !labels.insert(s.label).second)
    {
      throw std::runtime_error("Invalid/duplicate sample label");
    }
    samples.push_back(s);
  }
  if (samples.empty())
  {
    throw std::runtime_error("No samples");
  }
  const std::string refLabel = std::string(reference).empty() ? samples.front().label : reference;
  auto ri = std::find_if(samples.begin(), samples.end(),
                         [&](const Sample &s)
                         { return s.label == refLabel; });
  if (ri == samples.end())
  {
    throw std::runtime_error("Reference label not in manifest");
  }
  std::rotate(samples.begin(), ri, std::next(ri));
  if (gSystem->mkdir(output, true) != 0 && gSystem->AccessPathName(output))
  {
    throw std::runtime_error("Cannot create output directory");
  }
  const std::string out(output);
  // Explicitly refuse overwriting a previous analysis.
  if (!gSystem->AccessPathName((out + "/results.root").c_str()))
  {
    throw std::runtime_error("Output already exists: use a new directory");
  }
  std::unique_ptr<TFile> result(TFile::Open((out + "/results.root").c_str(), "CREATE"));
  if (!result || result->IsZombie())
  {
    throw std::runtime_error("Cannot create results.root");
  }
  std::ofstream available(out + "/availability.tsv");
  std::ofstream csv(out + "/slice_fits.csv");
  std::ofstream trends(out + "/trends.csv");
  available << "sample\thistogram\tfiles_present\tfiles_total\tstatus\n";
  csv << std::setprecision(12)
      << "sample,observable,bin,x,x_half_width,in_range_counts,status,root_status,covariance,mean,"
         "mean_error,"
         "sigma,sigma_error,fit_low,fit_high,chi2_ndf,underflow,overflow,ndf,background,fit_"
         "model\n";
  trends << std::setprecision(12) << "sample,observable,parameter,value,error\n";
  std::ofstream driftCSV(out + "/drift_summary.csv");
  std::ofstream summary(out + "/analysis_summary.tsv");
  summary << "sample\trdca_status\n";
  driftCSV << std::setprecision(12)
           << "sample,side,source,t0_valid,cm_valid,t0_bin,cm_bin,t0_width10_90,cm_width10_90,v_"
              "data_cm_us,v_"
              "phg_edges_cm_us,v_phg_baseline_cm_us,v_phg_adc_weighted_cm_us,z_pad_edge_cm,z_cm_"
              "edge_cm,z_"
              "pad_geometry_cm,cm_halfwidth_cm,time_bin_ns,raw_to_unpacked_shift_bins\n";
  const Config cfg{ptMin, ptMax, minEntries, phiGroup, etaGroup, prefix,
                   saveSlices, ratioFloor, phiQpt, rawTimeShift, fitBackground};
  folder(result.get(), "metadata")->cd();
  std::ostringstream config;
  config << "reference=" << refLabel << "; prefix=" << prefix << "; pt=[" << ptMin << ',' << ptMax
         << "); min_entries=" << minEntries << "; phi_group=" << phiGroup
         << "; eta_group=" << etaGroup
         << "; fits=two-pass Poisson/bin-integrated Gaussian; mean comparison=difference; other "
            "comparisons "
            "as named; errors=independent samples, fixed exposure";
  config << "; signed_ratio_floor_cm=" << ratioFloor
         << "; signed_ratio_reference_significance>3; save_slices=" << saveSlices
         << "; phi_qpt=" << phiQpt;
  config << "; fit_model="
         << (fitBackground ? "Gaussian+constant_background_modal_seed" : "Gaussian_quantile_seed");
  config << "; final_pass_quality_gate_only=true";
  config << "; raw_to_unpacked_shift=" << rawTimeShift << "; output_schema=2";
  TNamed metadata("configuration", config.str().c_str());
  metadata.Write();
  TNamed runConditions("run_conditions_json", conditions);
  runConditions.Write();
  Curves ref;
  Curves runTrends;
  std::vector<std::string> runLabels;
  for (const auto &sample : samples)
  {
    std::cout << "TPCDistortionQA: processing " << sample.label << std::endl;
    Hists hs = load(sample, cfg, available);
    DriftInput drift = TpcDriftQA::load(sample, available, cfg.rawTimeShift);
    TpcDriftQA::apply_reference(drift, sample, phgReferences);
    if (hs.empty() && !drift.time[0] && !drift.time[1])
    {
      throw std::runtime_error("No complete supported QA histograms for " + sample.label);
    }
    auto *dir = result->mkdir(sample.label.c_str());
    ObjectIndex outputIndex;
    for (const std::string name : {"rdca", "charge", "drift", "diagnostics", "metadata"})
    {
      folder(dir, name);
    }
    Curves curves = TpcDistortionQA::analyze(hs, sample, cfg, dir, csv, trends, outputIndex);
    bool hasDca = false;
    bool hasFit = false;
    for (const auto &entry : hs)
    {
      if (entry.first.find("dcaxyorigin_phi") != std::string::npos ||
          entry.first.find("_rdca_") != std::string::npos)
      {
        hasDca = true;
      }
    }
    for (const auto &entry : curves)
    {
      if (ends(entry.first, "_mean") && !entry.second.x.empty())
      {
        hasFit = true;
      }
    }
    std::string dcaState =
        "No supported DCA input histograms found; check histogram "
        "prefix and availability.tsv";
    if (hasDca)
    {
      dcaState = hasFit ? "DCA inputs found; Gaussian profiles available"
                        : "DCA inputs found, but no accepted Gaussian slices; "
                          "see diagnostics and min-entries";
    }
    std::cout << sample.label << ": " << dcaState << std::endl;
    summary << sample.label << '\t' << dcaState << '\n';
    summary.flush();
    folder(dir, "metadata")->cd();
    TNamed dcaAvailability("rdca_availability", dcaState.c_str());
    dcaAvailability.Write();
    Curves driftCurves = TpcDriftQA::analyze(drift, sample, dir, driftCSV, outputIndex);
    for (const auto &entry : driftCurves)
    {
      if (entry.second.xTitle != "side (0 South, 1 North)")
      {
        continue;
      }
      for (size_t i = 0; i < entry.second.x.size(); ++i)
      {
        const std::string key = entry.first + "_side" + std::to_string(int(entry.second.x[i]));
        auto &trend = runTrends[key];
        trend.title = key;
        trend.xTitle = "sample index (see bin labels)";
        trend.yTitle = entry.second.yTitle;
        point(trend, runLabels.size() + .5, .5, entry.second.y[i], entry.second.ey[i]);
      }
    }
    runLabels.push_back(sample.label);
    curves.insert(driftCurves.begin(), driftCurves.end());
    if (sample.label == refLabel)
    {
      ref = curves;
    }
    if (sample.label != refLabel)
    {
      for (const auto &entry : curves)
      {
        auto found = ref.find(entry.first);
        if (found == ref.end())
        {
          continue;
        }
        writeCurve(combine(entry.second, found->second, entry.second.comparison),
                   entry.first + "_vs_reference", dir, outputIndex);
        if (entry.first.find("_mean") != std::string::npos ||
            entry.first.find("_charge_even") != std::string::npos ||
            entry.first.find("_charge_odd") != std::string::npos)
        {
          writeCurve(combine(entry.second, found->second, "signed_ratio", cfg.ratioFloor),
                     entry.first + "_ratio_to_reference", dir, outputIndex);
        }
      }
    }
    if (makePlots)
    {
      plot(curves, ref, sample, refLabel, out + "/" + sample.label + ".pdf", cfg);
      plotDcaDistributions(dir, out + "/" + sample.label + "_rdca_distributions.pdf");
    }
    std::ostringstream index;
    index << "internal_key\tROOT_path_relative_to_run\n";
    for (const auto &entry : outputIndex)
    {
      index << entry.first << '\t' << entry.second << '\n';
    }
    folder(dir, "metadata")->cd();
    TNamed objectIndex("object_index", index.str().c_str());
    objectIndex.Write();
  }
  if (!runTrends.empty())
  {
    auto *dir = folder(result.get(), "run_trends/drift");
    TCanvas overview("drift_overview", "Drift QA run trends", 1400, 700);
    const std::string pdf = out + "/drift_run_trends.pdf";
    if (makePlots)
    {
      overview.Print((pdf + "[").c_str());
    }
    ObjectIndex trendIndex;
    for (auto &entry : runTrends)
    {
      auto &curve = entry.second;
      for (size_t i = 0; i <= runLabels.size(); ++i)
      {
        curve.edges.push_back(i);
      }
      const auto saved = writeCurve(curve, entry.first, dir, trendIndex, true);
      std::unique_ptr<TH1> h(dynamic_cast<TH1 *>(dir->Get(saved.c_str())));
      if (!h)
      {
        throw std::runtime_error("Missing run-trend histogram: " + saved);
      }
      h->SetDirectory(nullptr);
      for (size_t i = 0; i < runLabels.size(); ++i)
      {
        h->GetXaxis()->SetBinLabel(i + 1, runLabels[i].c_str());
      }
      dir->cd();
      h->Write("", TObject::kOverwrite);
      auto g = graph(curve, entry.first);
      g->SetMarkerStyle(20);
      g->Draw("AP");
      if (makePlots)
      {
        overview.Print(pdf.c_str());
      }
    }
    for (int side = 0; side < 2; ++side)
    {
      std::vector<std::unique_ptr<TGraphErrors>> graphs;
      TLegend legend(.60, .72, .98, .93);
      double low = std::numeric_limits<double>::infinity();
      double high = -low;
      const std::vector<std::string> names = {"data_velocity", "phg_edge_velocity",
                                              "phg_baseline_velocity", "phg_adc_weighted_velocity"};
      for (size_t index = 0; index < names.size(); ++index)
      {
        auto it = runTrends.find("drift_" + names[index] + "_side" + std::to_string(side));
        if (it == runTrends.end() || it->second.x.empty())
        {
          continue;
        }
        for (double y : it->second.y)
        {
          low = std::min(low, y);
          high = std::max(high, y);
        }
        auto g = graph(it->second, names[index]);
        const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
        g->SetMarkerColor(colors[index]);
        g->SetLineColor(colors[index]);
        g->SetMarkerStyle(20 + index);
        legend.AddEntry(g.get(), names[index].c_str(), "lp");
        graphs.push_back(std::move(g));
      }
      if (!graphs.empty())
      {
        overview.Clear();
        graphs.front()->SetTitle(("Data and PHGarfield velocities, side " + std::to_string(side) +
                                  ";sample index;v [cm/#mus]")
                                     .c_str());
        const double margin = std::max(.005, .2 * (high - low));
        graphs.front()->SetMinimum(low - margin);
        graphs.front()->SetMaximum(high + margin);
        graphs.front()->Draw("AP");
        for (size_t i = 1; i < graphs.size(); ++i)
        {
          graphs[i]->Draw("P SAME");
        }
        legend.Draw();
        if (makePlots)
        {
          overview.Print(pdf.c_str());
        }
      }
    }
    if (makePlots)
    {
      overview.Print((pdf + "]").c_str());
    }
  }
  result->Close();
  std::ofstream done(out + "/SUCCESS");
  done << "Analysis completed\n";
  std::cout << "TPCDistortionQA: wrote " << out << std::endl;
  return 0;
}
