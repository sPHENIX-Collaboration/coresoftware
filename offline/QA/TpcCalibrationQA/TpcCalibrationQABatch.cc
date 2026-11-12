#include "TpcCalibrationQA.h"
#include "TpcQAInternal.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2.h>
#include <TKey.h>
#include <TList.h>
#include <TNamed.h>
#include <TROOT.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Batch helpers: one pass through source segments, then portable summary
// collection.

namespace TpcCalibrationQABatch
{
  using TpcCalibrationQADetail::combine;
  using TpcCalibrationQADetail::Config;
  using TpcCalibrationQADetail::Curve;
  using TpcCalibrationQADetail::Curves;
  using TpcCalibrationQADetail::edges;
  using TpcCalibrationQADetail::ends;
  using TpcCalibrationQADetail::folder;
  using TpcCalibrationQADetail::hasKey;
  using TpcCalibrationQADetail::Hists;
  using TpcCalibrationQADetail::location;
  using TpcCalibrationQADetail::ObjectIndex;
  using TpcCalibrationQADetail::plain;
  using TpcCalibrationQADetail::plot;
  using TpcCalibrationQADetail::point;
  using TpcCalibrationQADetail::sameBins;
  using TpcCalibrationQADetail::Sample;
  using TpcCalibrationQADetail::wanted;
  using TpcCalibrationQADetail::writeCurve;
  bool keepInput(const std::string &name, const std::string &prefix)
  {
    for (const auto &key : wanted())
    {
      if (name == prefix + key)
      {
        return true;
      }
    }
    return name.rfind("h_TpcRawHitQA_nhits_sample_sec", 0) == 0 ||
           name.rfind("h_TpcDriftVelocityQA_", 0) == 0 ||
           name.rfind("h_PHGarfieldDriftVelocityQA_", 0) == 0;
  }
  void cache(const char *listing, const char *output, const char *prefix)
  {
    TH1::AddDirectory(false);
    std::ifstream input(listing);
    if (!input)
    {
      throw std::runtime_error("Cannot read cache input list");
    }
    Hists sums;
    std::set<std::string> schema;
    std::set<std::string> seen;
    std::string path;
    int count = 0;
    while (std::getline(input, path))
    {
      if (path.empty() || path[0] == '#')
      {
        continue;
      }
      if (!seen.insert(path).second)
      {
        throw std::runtime_error("Repeated cache input");
      }
      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
      {
        throw std::runtime_error("Bad cache input: " + path);
      }
      std::set<std::string> present;
      TIter next(file->GetListOfKeys());
      while (auto *key = dynamic_cast<TKey *>(next()))
      {
        std::string name = key->GetName();
        if (!keepInput(name, prefix) || hasKey(present, name))
        {
          continue;
        }
        auto *h = dynamic_cast<TH1 *>(file->Get(name.c_str()));
        if (!h)
        {
          continue;
        }
        present.insert(name);
        if (hasKey(sums, name))
        {
          if (!sameBins(sums[name].get(), h))
          {
            throw std::runtime_error("Mixed cache binning: " + name);
          }
          // Geometry and time conventions must stay constant when calibration
          // sums are merged.
          if (name == "h_TpcDriftVelocityQA_h_calibration_sums")
          {
            for (int bin : {6, 7, 8, 9, 10})
            {
              double n = sums[name]->GetBinContent(1);
              double m = h->GetBinContent(1);
              if (!(n > 0 && m > 0) ||
                  std::abs(sums[name]->GetBinContent(bin) / n - h->GetBinContent(bin) / m) > 1e-6)
              {
                throw std::runtime_error("Mixed drift geometry/time settings in cache");
              }
            }
          }
          if (!sums[name]->Add(h))
          {
            throw std::runtime_error("Failed histogram addition");
          }
        }
        else
        {
          sums[name].reset(dynamic_cast<TH1 *>(h->Clone(name.c_str())));
          sums[name]->SetDirectory(nullptr);
        }
      }
      if (count && schema != present)
      {
        throw std::runtime_error(
            "Mixed supported QA schemas within a run; split "
            "the inventory by reconstruction configuration");
      }
      schema = present;
      ++count;
      if (count % 1000 == 0)
      {
        std::cout << "Cached " << count << " segments" << std::endl;
      }
    }
    if (!count || sums.empty())
    {
      throw std::runtime_error("No supported histograms to cache");
    }
    TFile out(output, "RECREATE");
    if (out.IsZombie())
    {
      throw std::runtime_error("Cannot create QA cache");
    }
    for (auto &entry : sums)
    {
      entry.second->Write();
    }
    TNamed note("cache_provenance",
                ("Selected QA only; source segments=" + std::to_string(count)).c_str());
    note.Write();
    out.Close();
  }

  Curves readCurves(TDirectory *dir)
  {
    auto *index = dynamic_cast<TNamed *>(dir->Get("metadata/object_index"));
    if (!index)
    {
      throw std::runtime_error("Missing summary object index");
    }
    Curves curves;
    std::istringstream lines(index->GetTitle());
    std::string line;
    std::getline(lines, line);
    while (std::getline(lines, line))
    {
      auto tab = line.find('\t');
      if (tab == std::string::npos)
      {
        continue;
      }
      std::string key = line.substr(0, tab);
      std::string path = line.substr(tab + 1);
      if (key.find("_vs_reference") != std::string::npos ||
          key.find("_ratio_to_reference") != std::string::npos)
      {
        continue;
      }
      auto *object = dir->Get(path.c_str());
      auto *mode = dynamic_cast<TNamed *>(dir->Get(("metadata/comparisons/" + key).c_str()));
      if (!mode)
      {
        continue;  // index also contains input TH2 maps and fit-result objects
      }
      Curve c;
      c.title = object ? object->GetTitle() : key;
      c.comparison = mode->GetTitle();
      if (auto *h = dynamic_cast<TH1 *>(object))
      {
        if (h->GetDimension() != 1)
        {
          continue;
        }
        auto *valid = dynamic_cast<TH1 *>(dir->Get(("diagnostics/" + path + "_valid").c_str()));
        if (path.rfind("diagnostics/", 0) == 0)
        {
          valid = dynamic_cast<TH1 *>(dir->Get((path + "_valid").c_str()));
        }
        if (!valid)
        {
          throw std::runtime_error("Missing validity mask: " + path);
        }
        c.xTitle = h->GetXaxis()->GetTitle();
        c.yTitle = h->GetYaxis()->GetTitle();
        c.edges = edges(h->GetXaxis());
        for (int b = 1; b <= h->GetNbinsX(); ++b)
        {
          if (valid->GetBinContent(b) > 0)
          {
            point(c, h->GetBinCenter(b), .5 * h->GetBinWidth(b), h->GetBinContent(b),
                  h->GetBinError(b));
          }
        }
      }
      else if (auto *g = dynamic_cast<TGraphErrors *>(object))
      {
        c.xTitle = g->GetXaxis()->GetTitle();
        c.yTitle = g->GetYaxis()->GetTitle();
        for (int i = 0; i < g->GetN(); ++i)
        {
          point(c, g->GetPointX(i), g->GetErrorX(i), g->GetPointY(i), g->GetErrorY(i));
        }
      }
      else
      {
        continue;
      }
      curves[key] = std::move(c);
    }
    return curves;
  }

  void copy(TDirectory *src, TDirectory *dst, bool maps)
  {
    // Keep a frame per directory to preserve depth-first key order without recursion.
    struct Frame
    {
      TDirectory *destination;
      TIter next;
      std::set<std::string> seen;
    };
    std::vector<Frame> stack;
    stack.push_back({dst, TIter(src->GetListOfKeys()), {}});
    while (!stack.empty())
    {
      auto &frame = stack.back();
      auto *key = dynamic_cast<TKey *>(frame.next());
      if (!key)
      {
        stack.pop_back();
        continue;
      }
      const std::string name = key->GetName();
      if (!frame.seen.insert(name).second)
      {
        continue;
      }
      // Read from the directory owning this key, as in the recursive traversal.
      auto *object = key->GetMotherDir()->Get(name.c_str());
      if (auto *sub = dynamic_cast<TDirectory *>(object))
      {
        auto *destination = folder(frame.destination, name);
        stack.push_back({destination, TIter(sub->GetListOfKeys()), {}});
        continue;
      }
      if (!maps)
      {
        if (auto *h = dynamic_cast<TH1 *>(object))
        {
          if (h->GetDimension() > 1)
          {
            continue;
          }
        }
        if (ends(name, "_valid_points"))
        {
          continue;
        }
      }
      frame.destination->cd();
      object->Write(name.c_str());
    }
  }

  struct Item
  {
    std::string group, window, label, path;
  };
  void collect(const char *manifest, const char *references, const char *output, bool maps = false,
               bool plots = false, double ratioFloor = .02)
  {
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    std::ifstream in(manifest);
    std::ifstream refs(references);
    if (!in || !refs)
    {
      throw std::runtime_error("Cannot read collection manifest/reference list");
    }
    std::vector<Item> items;
    std::string line;
    while (std::getline(in, line))
    {
      std::istringstream row(line);
      Item i;
      std::getline(row, i.group, '\t');
      std::getline(row, i.window, '\t');
      std::getline(row, i.label, '\t');
      std::getline(row, i.path);
      items.push_back(i);
    }
    std::map<std::string, std::string> reference;
    while (std::getline(refs, line))
    {
      auto p = line.find('\t');
      reference[line.substr(0, p)] = line.substr(p + 1);
    }
    std::map<std::string, Curves> refCurves;
    for (const auto &i : items)
    {
      if (reference[i.group] == i.label)
      {
        TFile file(i.path.c_str(), "READ");
        auto *dir = file.GetDirectory(i.label.c_str());
        if (file.IsZombie() || !dir)
        {
          throw std::runtime_error("Bad reference summary");
        }
        refCurves[i.group + "/" + i.window] = readCurves(dir);
      }
    }
    TFile result(output, "RECREATE");
    if (result.IsZombie())
    {
      throw std::runtime_error("Cannot create combined summary");
    }
    std::map<std::string, std::map<std::string, std::pair<double, double>>> trends;
    for (const auto &i : items)
    {
      const std::string base = i.group + "/" + i.window;
      if (!hasKey(refCurves, base))
      {
        throw std::runtime_error("Reference missing for " + base);
      }
      TFile file(i.path.c_str(), "READ");
      auto *source = file.GetDirectory(i.label.c_str());
      if (file.IsZombie() || !source)
      {
        throw std::runtime_error("Bad per-run summary");
      }
      // Accepted harmonic/linear fit parameters become small run-trend histograms
      // too.
      std::function<void(TDirectory *, const std::string &)> extractFits;
      extractFits = [&](TDirectory *current, const std::string &path)
      {
        TIter next(current->GetListOfKeys());
        std::set<std::string> seen;
        while (auto *key = dynamic_cast<TKey *>(next()))
        {
          std::string name = key->GetName();
          if (!seen.insert(name).second)
          {
            continue;
          }
          auto *object = current->Get(name.c_str());
          if (auto *sub = dynamic_cast<TDirectory *>(object))
          {
            extractFits(sub, path + name + "/");
          }
          else if (auto *fit = dynamic_cast<TFitResult *>(object))
          {
            if (fit->Status() != 0 || fit->CovMatrixStatus() < 2)
            {
              continue;
            }
            for (unsigned p = 0; p < fit->NPar(); ++p)
            {
              std::string title = path + name + "_parameter_" + std::to_string(p);
              std::replace(title.begin(), title.end(), '/', '_');
              trends[std::string(base).append("/run_trends/rdca/").append(title)][i.label] = {
                  fit->Parameter(p), fit->ParError(p)};
            }
          }
        }
      };
      if (auto *diagnostics = source->GetDirectory("diagnostics"))
      {
        extractFits(diagnostics, "");
      }
      auto curves = readCurves(source);
      auto *destination = folder(&result, base + "/runs/" + i.label);
      copy(source, destination, maps);
      auto *worker = folder(destination, "metadata/analysis");
      if (auto *global = file.GetDirectory("metadata"))
      {
        copy(global, worker, true);
      }
      ObjectIndex outputIndex;
      for (const auto &entry : curves)
      {
        const auto &c = entry.second;
        if (entry.first.rfind("drift_", 0) == 0 && c.xTitle == "side (0 South, 1 North)")
        {
          for (size_t j = 0; j < c.x.size(); ++j)
          {
            trends[base + "/run_trends/drift/" + location(entry.first).name + "_" +
                   (int(c.x[j]) == 0 ? "south" : "north")][i.label] = {c.y[j], c.ey[j]};
          }
        }
        auto found = refCurves[base].find(entry.first);
        if (i.label == reference[i.group] || found == refCurves[base].end())
        {
          continue;
        }
        writeCurve(combine(c, found->second, c.comparison), entry.first + "_vs_reference",
                   destination, outputIndex);
        if (entry.first.find("_mean") != std::string::npos ||
            entry.first.find("_charge_even") != std::string::npos ||
            entry.first.find("_charge_odd") != std::string::npos)
        {
          writeCurve(combine(c, found->second, "signed_ratio", ratioFloor),
                     entry.first + "_ratio_to_reference", destination, outputIndex);
        }
      }
      // Keep an explicit index of comparisons; original index describes source
      // observables.
      folder(destination, "metadata")->cd();
      std::ostringstream index;
      for (const auto &pair : outputIndex)
      {
        index << pair.first << '\t' << pair.second << '\n';
      }
      TNamed comparisonIndex("reference_comparison_index", index.str().c_str());
      comparisonIndex.Write();
      TNamed referenceNote("reference_run", reference[i.group].c_str());
      referenceNote.Write();
      if (plots)
      {
        Config cfg{.5, 1., 100, 3, 4, "h_TpcSeedsQA_", false, ratioFloor, false, 40.};
        Sample sample{i.label, "", ""};
        std::string pdf =
            std::string(output) + "_" + i.group + "_" + i.window + "_" + i.label + ".pdf";
        std::replace(pdf.begin() + std::string(output).size(), pdf.end(), '/', '_');
        plot(curves, refCurves[base], sample, reference[i.group], pdf, cfg);
      }
    }
    for (const auto &entry : trends)
    {
      auto p = entry.first.rfind('/');
      auto *dir = folder(&result, entry.first.substr(0, p));
      std::string name = entry.first.substr(p + 1);
      TH1D h(name.c_str(), (plain(name) + ";run;value (errors as provided)").c_str(),
             entry.second.size(), 0, entry.second.size());
      h.SetDirectory(nullptr);
      h.SetStats(false);
      h.SetOption("E1");
      int bin = 0;
      for (const auto &v : entry.second)
      {
        ++bin;
        h.GetXaxis()->SetBinLabel(bin, v.first.c_str());
        h.SetBinContent(bin, v.second.first);
        h.SetBinError(bin, v.second.second);
      }
      dir->cd();
      h.Write();
    }
    folder(&result, "metadata")->cd();
    TNamed scope("summary_scope",
                 maps ? "Includes input maps"
                      : "Compact: 1D profiles, validity masks, fit results, metadata; "
                        "2D/3D input maps remain in per-run SDCC outputs");
    scope.Write();
    TNamed limits("interpretation",
                  "Observable/reference ratios are not k_eff ratios. Derive "
                  "response using controlled k_eff variations; preserve side, "
                  "charge, selection and reconstruction configuration. Drift "
                  "edge uncertainties are not estimated.");
    limits.Write();
    result.Close();
  }

  // Attach readable side scalars and run trends to the portable ROOT summary.
  void attachConditions(const char *table, const char *output)
  {
    using TpcCalibrationQADetail::folder;
    std::ifstream input(table);
    if (!input)
    {
      throw std::runtime_error("Cannot read condition scalar table");
    }
    std::map<std::string, std::map<int, std::vector<double>>> rows;
    std::string line;
    while (std::getline(input, line))
    {
      std::istringstream fields(line);
      std::vector<std::string> tokens;
      std::string field;
      while (std::getline(fields, field, '\t'))
      {
        tokens.push_back(field);
      }
      if (tokens.size() != 7)
      {
        throw std::runtime_error("Invalid condition scalar row");
      }
      const int side = std::stoi(tokens[3]);
      if (side < 0 || side > 1)
      {
        throw std::runtime_error("Invalid condition scalar side");
      }
      const std::string key = tokens[0] + "/" + tokens[1] + "/runs/" + tokens[2];
      if (hasKey(rows[key], side))
      {
        throw std::runtime_error("Duplicate condition scalar side");
      }
      rows[key][side] = {std::stod(tokens[4]), std::stod(tokens[5]), std::stod(tokens[6])};
    }
    TFile file(output, "UPDATE");
    if (file.IsZombie())
    {
      throw std::runtime_error("Cannot update condition scalars");
    }
    const std::vector<std::string> names = {"gem_r1_current_run_average", "keff_base_snapshot",
                                            "keff_applied"};
    const std::vector<std::string> units = {"R1 current [#muA]", "base k_{eff} snapshot",
                                            "verified applied k_{eff}"};
    std::map<std::string, std::map<std::string, double>> trends;
    for (const auto &run : rows)
    {
      auto cut = run.first.rfind("/runs/");
      const auto group = run.first.substr(0, cut);
      const auto label = run.first.substr(cut + 6);
      auto *dir = folder(&file, run.first + "/conditions");
      for (std::size_t metric = 0; metric < names.size(); ++metric)
      {
        bool any = false;
        for (const auto &side : run.second)
        {
          any |= std::isfinite(side.second[metric]);
        }
        if (!any)
        {
          continue;
        }
        const std::string name = names[metric] + "_vs_side";
        TH1D h(name.c_str(), (names[metric] + ";side (0 South, 1 North);" + units[metric]).c_str(), 2,
               -.5, 1.5);
        TH1D mask((name + "_valid").c_str(), "Validity: 1 means recorded;side;valid", 2, -.5, 1.5);
        h.SetDirectory(nullptr);
        mask.SetDirectory(nullptr);
        h.SetStats(false);
        h.SetOption("P");
        for (const auto &side : run.second)
        {
          if (std::isfinite(side.second[metric]))
          {
            h.SetBinContent(side.first + 1, side.second[metric]);
            mask.SetBinContent(side.first + 1, 1);
            trends[group + "/run_trends/conditions/" + names[metric] +
                   (side.first == 0 ? "_south" : "_north")][label] = side.second[metric];
          }
        }
        dir->cd();
        h.Write("", TObject::kOverwrite);
        folder(&file, run.first + "/diagnostics/conditions")->cd();
        mask.Write("", TObject::kOverwrite);
      }
      folder(&file, run.first + "/metadata")->cd();
      TNamed scope("conditions_scope",
                   "See analysis/run_conditions_json for payload, side, units, coverage and kind. "
                   "CDB base snapshots and derived fallback coefficients are not verified historical "
                   "applied k_eff. Missing applied values have no histogram.");
      scope.Write("", TObject::kOverwrite);
    }
    for (const auto &entry : trends)
    {
      const auto slash = entry.first.rfind('/');
      const auto name = entry.first.substr(slash + 1);
      TH1D h(name.c_str(), (name + ";run;recorded condition").c_str(), entry.second.size(), 0,
             entry.second.size());
      h.SetDirectory(nullptr);
      h.SetStats(false);
      h.SetOption("P");
      int bin = 0;
      for (const auto &value : entry.second)
      {
        ++bin;
        h.GetXaxis()->SetBinLabel(bin, value.first.c_str());
        h.SetBinContent(bin, value.second);
      }
      folder(&file, entry.first.substr(0, slash))->cd();
      h.Write("", TObject::kOverwrite);
    }
    file.Close();
  }
}  // namespace TpcCalibrationQABatch

void TpcCalibrationQA::cache(const char *listing, const char *output, const char *prefix)
{
  TpcCalibrationQABatch::cache(listing, output, prefix);
}
void TpcCalibrationQA::collect(const char *manifest, const char *references, const char *output,
                               bool maps, bool plots, double ratioFloor)
{
  TpcCalibrationQABatch::collect(manifest, references, output, maps, plots, ratioFloor);
}
void TpcCalibrationQA::attach_conditions(const char *table, const char *output)
{
  TpcCalibrationQABatch::attachConditions(table, output);
}
