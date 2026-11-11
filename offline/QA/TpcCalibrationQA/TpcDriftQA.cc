#include "TpcDriftQA.h"

namespace TpcCalibrationQADetail
{
  void TpcDriftQA::apply_reference(DriftInput &d, const Sample &sample, const char *path)
  {
    if (!path || !*path)
    {
      return;
    }
    std::ifstream input(path);
    if (!input)
    {
      throw std::runtime_error("Cannot open PHGarfield reference table");
    }
    std::string line;
    while (std::getline(input, line))
    {
      std::istringstream row(line);
      std::string label;
      std::string side;
      std::string velocity;
      std::string t0;
      std::string source;
      std::getline(row, label, '\t');
      std::getline(row, side, '\t');
      std::getline(row, velocity, '\t');
      std::getline(row, t0, '\t');
      std::getline(row, source, '\t');
      std::string clock;
      std::string pad;
      std::string half;
      std::getline(row, clock, '\t');
      std::getline(row, pad, '\t');
      std::getline(row, half);
      if (label == sample.label)
      {
        for (const auto &field : {std::make_pair(clock, &d.clock), std::make_pair(pad, &d.pad), std::make_pair(half, &d.cmHalf)})
        {
          if (field.first.empty())
          {
            continue;
          }
          double value = std::stod(field.first);
          if (!std::isfinite(value) || value <= 0)
          {
            throw std::runtime_error("Invalid external PHGarfield geometry/clock");
          }
          if (d.conditions && std::abs(*field.second - value) > 1e-6)
          {
            throw std::runtime_error("Native and external PHGarfield geometry/clock differ");
          }
          *field.second = value;
        }
      }
      if (label != sample.label)
      {
        continue;
      }
      int index = std::stoi(side);
      if (index < 0 || index > 1 || !d.referenceSource[index].empty())
      {
        throw std::runtime_error("Invalid/duplicate PHGarfield reference side");
      }
      d.referenceSource[index] = source;
      if (!velocity.empty())
      {
        d.referenceV[index] = std::stod(velocity);
      }
      if (!t0.empty())
      {
        d.referenceT0[index] = std::stod(t0);
      }
    }
  }

  double median(std::vector<double> x)
  {
    if (x.empty())
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(x.begin(), x.end());
    return x.size() % 2 ? x[x.size() / 2] : .5 * (x[x.size() / 2 - 1] + x[x.size() / 2]);
  }

  double windowMedian(const std::vector<double> &x, const std::vector<double> &y, double lo, double hi)
  {
    std::vector<double> values;
    for (size_t i = 0; i < x.size(); ++i)
    {
      if (x[i] >= lo && x[i] <= hi && std::isfinite(y[i]))
      {
        values.push_back(y[i]);
      }
    }
    return median(values);
  }

  Edge edge(const std::vector<double> &x, const std::vector<double> &raw, double ll, double lh, double rl,
            double rh, double sl, double sh, double expected, int direction, bool clean)
  {
    Edge e;
    if (x.size() != raw.size() || x.size() < 5)
    {
      return e;
    }
    auto y = raw;
    if (clean)
    {
      for (size_t i = 1; i + 1 < y.size(); ++i)
      {
        y[i] = median({raw[i - 1], raw[i], raw[i + 1]});
      }
    }
    e.left = windowMedian(x, y, ll, lh);
    e.right = windowMedian(x, y, rl, rh);
    const double height = e.right - e.left;
    if (!std::isfinite(height) || direction * height <= 0 ||
        std::abs(height) < .01 * std::max(std::abs(e.left), std::abs(e.right)))
    {
      return e;
    }
    std::vector<double> crossings;
    for (double fraction : {.1, .5, .9})
    {
      const double target = e.left + fraction * height;
      std::vector<double> found;
      for (size_t i = 0; i + 1 < x.size(); ++i)
      {
        if (x[i] < sl || x[i + 1] > sh || direction * (y[i + 1] - y[i]) <= 0)
        {
          continue;
        }
        if ((y[i] - target) * (y[i + 1] - target) <= 0)
        {
          found.push_back(x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]));
        }
      }
      // Unlike the old notebook, never turn a nearest-bin fallback into an edge.
      if (found.empty())
      {
        return e;
      }
      crossings.push_back(*std::min_element(found.begin(), found.end(), [&](double a, double b)
                                            { return std::abs(a - expected) < std::abs(b - expected); }));
    }
    e.x10 = crossings[0];
    e.x50 = crossings[1];
    e.x90 = crossings[2];
    e.valid = e.x10 <= e.x50 && e.x50 <= e.x90;
    return e;
  }

  Hist clone1(TH1 *h, const std::string &name)
  {
    if (!h || h->GetDimension() != 1)
    {
      return {};
    }
    Hist out(static_cast<TH1 *>(h->Clone(name.c_str())));
    out->SetDirectory(nullptr);
    return out;
  }

  Hist readTime(TFile *file, int side, std::string &source)
  {
    const std::string suffix = "_side" + std::to_string(side);
    if (auto h = clone1(dynamic_cast<TH1 *>(file->Get(("h_TpcDriftVelocityQA_h_adc_time" + suffix).c_str())),
                        "drift_time"))
    {
      source = "cleaned_ADC_side";
      return h;
    }
    if (auto *h =
            dynamic_cast<TH2 *>(file->Get(("h_PHGarfieldDriftVelocityQA_h_adc_tbin_vs_fee" + suffix).c_str())))
    {
      Hist sum(h->ProjectionY("drift_time", 1, h->GetNbinsX(), "e"));
      sum->SetDirectory(nullptr);
      source = "standalone_cleaned_ADC_FEE_sum";
      // These histograms were filled at integer raw bins into [0,1000): correct
      // their half-bin display offset. The CM-t0 interval is unchanged.
      if (sum->GetXaxis()->GetXmin() == 0 && sum->GetXaxis()->GetBinWidth(1) == 1)
      {
        sum->GetXaxis()->SetLimits(-.5, sum->GetNbinsX() - .5);
        sum->ResetStats();
      }
      return sum;
    }
    Hist total;
    // Raw hardware sectors 0..11 North, 12..23 South; thresholded sample counts.
    for (int sector = side == 0 ? 12 : 0; sector < (side == 0 ? 24 : 12); ++sector)
    {
      for (int region = 0; region < 3; ++region)
      {
        const std::string name =
            "h_TpcRawHitQA_nhits_sample_sec" + std::to_string(sector) + "_R" + std::to_string(region);
        auto *h = dynamic_cast<TH1 *>(file->Get(name.c_str()));
        if (!h || h->GetDimension() != 1)
        {
          return {};
        }
        if (!total)
        {
          total = clone1(h, "drift_time");
        }
        else
        {
          if (!sameBins(total.get(), h))
          {
            throw std::runtime_error("Raw timing bins differ");
          }
          total->Add(h);
        }
      }
    }
    if (total)
    {
      source = "production_thresholded_sample_counts";
    }
    return total;
  }

  Hist readZ(TFile *file, int side)
  {
    const std::string suffix = "_side" + std::to_string(side);
    if (auto h = clone1(dynamic_cast<TH1 *>(file->Get(("h_TpcDriftVelocityQA_h_adc_z" + suffix).c_str())),
                        "drift_z"))
    {
      return h;
    }
    if (auto *h =
            dynamic_cast<TH2 *>(file->Get(("h_PHGarfieldDriftVelocityQA_h_adc_z_vs_fee" + suffix).c_str())))
    {
      Hist sum(h->ProjectionY("drift_z", 1, h->GetNbinsX(), "e"));
      sum->SetDirectory(nullptr);
      return sum;
    }
    return {};
  }

  DriftInput TpcDriftQA::load(const Sample &sample, std::ostream &availability, double defaultRawShift)
  {
    DriftInput d;
    const std::string listing = sample.driftList.empty() ? sample.list : sample.driftList;
    std::ifstream in(listing);
    if (!in)
    {
      throw std::runtime_error("Cannot read drift list: " + listing);
    }
    int files = 0;
    int nt[2]{0, 0};
    int nz[2]{0, 0};
    int nc[2]{0, 0};
    int nw[2]{0, 0};
    std::set<std::string> seen;
    std::string path;
    bool geometrySeen = false;
    while (std::getline(in, path))
    {
      if (path.empty() || path[0] == '#')
      {
        continue;
      }
      if (!seen.insert(path).second)
      {
        throw std::runtime_error("Duplicate drift input");
      }
      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      if (!file || file->IsZombie() || file->TestBit(TFile::kRecovered))
      {
        throw std::runtime_error("Bad drift ROOT input");
      }
      ++files;
      if (auto *h = dynamic_cast<TH1 *>(file->Get("h_TpcDriftVelocityQA_h_calibration_sums")))
      {
        if (!d.conditions)
        {
          d.conditions = clone1(h, "drift_calibration_sums");
        }
        else
        {
          if (!sameBins(d.conditions.get(), h))
          {
            throw std::runtime_error("Calibration schemas differ");
          }
          d.conditions->Add(h);
        }
        const double count = h->GetBinContent(1);
        if (!(count > 0))
        {
          throw std::runtime_error("Invalid drift calibration count");
        }
        const double pad = h->GetBinContent(6) / count;
        const double clock = h->GetBinContent(7) / count;
        const double half = h->GetBinContent(10) / count;
        if (geometrySeen && (std::abs(d.pad - pad) > 1e-6 || std::abs(d.clock - clock) > 1e-6 ||
                             std::abs(d.cmHalf - half) > 1e-6))
        {
          throw std::runtime_error("Mixed drift geometry/clock within sample");
        }
        d.pad = pad;
        d.clock = clock;
        d.cmHalf = half;
        geometrySeen = true;
      }
      for (int side = 0; side < 2; ++side)
      {
        std::string source;
        Hist t = readTime(file.get(), side, source);
        Hist z = readZ(file.get(), side);
        if (t)
        {
          if (source == "production_thresholded_sample_counts")
          {
            const double shift = std::isfinite(sample.rawTimeShift) ? sample.rawTimeShift : defaultRawShift;
            if (!std::isfinite(shift))
            {
              throw std::runtime_error("Raw time shift must be finite");
            }
            auto raw = clone1(t.get(), "time_raw");
            if (!d.rawTime[side])
            {
              d.rawTime[side] = std::move(raw);
            }
            else
            {
              d.rawTime[side]->Add(raw.get());
            }
            d.appliedShift[side] = shift;
            const double low = t->GetXaxis()->GetXmin();
            const double high = t->GetXaxis()->GetXmax();
            t->GetXaxis()->SetLimits(low - shift, high - shift);
            t->ResetStats();
          }
          if (!d.source[side].empty() && d.source[side] != source)
          {
            throw std::runtime_error("Mixed timing selections within sample");
          }
          d.source[side] = source;
          ++nt[side];
          if (!d.time[side])
          {
            d.time[side] = std::move(t);
          }
          else
          {
            if (!sameBins(d.time[side].get(), t.get()))
            {
              throw std::runtime_error("Timing segment binning differs");
            }
            d.time[side]->Add(t.get());
          }
        }
        if (z)
        {
          ++nz[side];
          if (!d.z[side])
          {
            d.z[side] = std::move(z);
          }
          else
          {
            if (!sameBins(d.z[side].get(), z.get()))
            {
              throw std::runtime_error("Z segment binning differs");
            }
            d.z[side]->Add(z.get());
          }
        }
        const std::string suffix = "_side" + std::to_string(side);
        auto *v = dynamic_cast<TH1 *>(file->Get(("h_TpcDriftVelocityQA_h_phg_mean_velocity" + suffix).c_str()));
        auto *c = dynamic_cast<TH1 *>(file->Get("h_TpcDriftVelocityQA_h_calibration_sums"));
        if (v && c)
        {
          ++nc[side];
          d.phgMean[side] += v->GetBinContent(1);
          d.calibrations[side] += c->GetBinContent(1);
        }
        else
        {
          auto *old =
              dynamic_cast<TH1 *>(file->Get(("h_PHGarfieldDriftVelocityQA_h_phg_vdrift" + suffix).c_str()));
          auto *n = dynamic_cast<TH1 *>(file->Get("h_PHGarfieldDriftVelocityQA_h_calibration_count"));
          if (old && n && n->GetBinContent(1) > 0)
          {
            ++nc[side];
            d.phgMean[side] += old->Integral() / old->GetNbinsX();
            d.calibrations[side] += n->GetBinContent(1);
          }
        }
        if (auto *w =
                dynamic_cast<TH1 *>(file->Get(("h_TpcDriftVelocityQA_h_adc_velocity_sums" + suffix).c_str())))
        {
          ++nw[side];
          d.weightedV[side] += w->GetBinContent(1);
          d.weightedADC[side] += w->GetBinContent(2);
        }
      }
    }
    for (int side = 0; side < 2; ++side)
    {
      for (const std::string kind : {"time", "phg_z"})
      {
        int count = kind == "time" ? nt[side] : nz[side];
        std::string state = "absent";
        if (files && count == files)
        {
          state = "complete";
        }
        else if (count)
        {
          state = "partial_SKIPPED";
        }
        availability << sample.label << "\tdrift_" << kind << "_side" << side << '\t' << count
                     << '\t' << files << '\t' << state << '\n';
      }
      if (!files || nt[side] != files)
      {
        d.time[side].reset();
        d.rawTime[side].reset();
      }
      if (!files || nz[side] != files)
      {
        d.z[side].reset();
      }
      if (nc[side] != files)
      {
        d.phgMean[side] = 0.;
        d.calibrations[side] = 0.;
      }
      if (nw[side] != files)
      {
        d.weightedV[side] = 0.;
        d.weightedADC[side] = 0.;
      }
    }
    return d;
  }

  std::pair<std::vector<double>, std::vector<double>> xy(const TH1 *h, int side, bool fold)
  {
    std::vector<std::pair<double, double>> points;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
      points.emplace_back((fold && side == 0 ? -1. : 1.) * h->GetBinCenter(i), h->GetBinContent(i));
    }
    std::sort(points.begin(), points.end());
    std::vector<double> x;
    std::vector<double> y;
    for (const auto &p : points)
    {
      x.push_back(p.first);
      y.push_back(p.second);
    }
    return {x, y};
  }

  std::pair<std::vector<double>, std::vector<double>> physicalZ(const TH1 *h, int side, double pad)
  {
    const auto points = xy(h, side, true);
    const auto &x = points.first;
    const auto &y = points.second;
    std::vector<double> xx;
    std::vector<double> yy;
    std::vector<size_t> regions;
    // Match notebook upper-branch selection in ~1 cm windows. Never bridge CM
    // or the pad plane, and do not average the downward binning aliases.
    const std::vector<double> limits = {x.front() - 1., 0., pad, x.back() + 1.};
    for (size_t region = 0; region + 1 < limits.size(); ++region)
    {
      std::vector<size_t> ids;
      for (size_t i = 0; i < x.size(); ++i)
      {
        if (x[i] >= limits[region] && x[i] < limits[region + 1])
        {
          ids.push_back(i);
        }
      }
      if (ids.empty())
      {
        continue;
      }
      for (double lo = x[ids.front()]; lo <= x[ids.back()]; lo += 1.05)
      {
        size_t best = x.size();
        for (size_t i : ids)
        {
          if (x[i] >= lo && x[i] < lo + 1.05 && (best == x.size() || y[i] > y[best]))
          {
            best = i;
          }
        }
        if (best < x.size())
        {
          xx.push_back(x[best]);
          yy.push_back(y[best]);
          regions.push_back(region);
        }
      }
    }
    auto cleanX = xx;
    auto cleanY = yy;
    cleanX.clear();
    cleanY.clear();
    for (size_t i = 0; i < xx.size(); ++i)
    {
      bool keep = true;
      if (i >= 2 && i + 2 < xx.size() && regions[i - 2] == regions[i] && regions[i + 2] == regions[i])
      {
        const double ref = median({yy[i - 2], yy[i - 1], yy[i + 1], yy[i + 2]});
        if (ref > 0 && yy[i] < .8 * ref)
        {
          keep = false;
        }
      }
      if (keep)
      {
        cleanX.push_back(xx[i]);
        cleanY.push_back(yy[i]);
      }
    }
    return {cleanX, cleanY};
  }

  Curves TpcDriftQA::analyze(DriftInput &d, const Sample &sample, TDirectory *dir, std::ostream &csv)
  {
    Curves curves;
    folder(dir, "metadata")->cd();
    if (d.conditions)
    {
      d.conditions->Write("drift_calibration_sums");
    }
    std::ostringstream geometryDescription;
    geometryDescription << "Timing geometry used: z_pad=" << d.pad << " cm; CM_half=" << d.cmHalf
                        << " cm; dt=" << d.clock << " ns. Native calibration sums take precedence; otherwise "
                                                    "supplied reference geometry/clock, then legacy defaults. Timing selection and "
                                                    "external reference records are stored separately.";
    TNamed convention("drift_convention", geometryDescription.str().c_str());
    convention.Write();
    const auto profile = [&](const std::string &name, int side, const std::vector<double> &x,
                             const std::vector<double> &y, const Edge &e, double lo, double hi,
                             const std::string &xTitle)
    {
      if (!e.valid || e.right == e.left)
      {
        return;
      }
      Curve c = makeCurve(sample.label + ": " + name + ", side " + std::to_string(side) + ", " + d.source[side], xTitle,
                          "normalized local edge");
      for (size_t i = 0; i < x.size(); ++i)
      {
        if (x[i] >= lo && x[i] <= hi)
        {
          point(c, x[i], 0., (y[i] - e.left) / (e.right - e.left), 0.);
        }
      }
      curves["drift_" + name + "_profile_side" + std::to_string(side)] = std::move(c);
    };
    const auto metric = [&](const std::string &name, int side, double value, const std::string &unit,
                            const std::string &comparison = "difference")
    {
      auto &c = curves["drift_" + name];
      c.title = sample.label + ": side-integrated drift QA (midpoint edges, errors not estimated)";
      c.xTitle = "side (0 South, 1 North)";
      c.yTitle = unit;
      c.comparison = comparison;
      c.edges = {-.5, .5, 1.5};
      point(c, side, .5, value, 0.);
    };
    for (int side = 0; side < 2; ++side)
    {
      const bool native = d.calibrations[side] > 0;
      folder(dir, "metadata")->cd();
      const std::string note = "native velocity: " + std::string(native ? "available" : "missing") +
                               "; external reference: " + (d.referenceSource[side].empty() ? "missing" : "provided; see phgarfield_reference_record_side" + std::to_string(side));
      TNamed availability(("phgarfield_availability_side" + std::to_string(side)).c_str(), note.c_str());
      availability.Write();
      if (!d.referenceSource[side].empty())
      {
        TNamed record(("phgarfield_reference_record_side" + std::to_string(side)).c_str(), d.referenceSource[side].c_str());
        record.Write();
      }
      if (native)
      {
        metric("phg_baseline_velocity", side, d.phgMean[side] / d.calibrations[side],
               "equal-FEE mean PHGarfield baseline [cm/#mus]", "ratio");
      }
      if (std::isfinite(d.referenceV[side]))
      {
        metric("phg_reference_velocity", side, d.referenceV[side],
               "supplied PHGarfield reference v [cm/#mus]", "ratio");
      }
      double configuredT0 = d.referenceT0[side];
      if (!std::isfinite(configuredT0) && d.conditions && d.conditions->GetBinContent(1) > 0)
      {
        configuredT0 = d.conditions->GetBinContent(8) / d.conditions->GetBinContent(1);
      }
      if (std::isfinite(configuredT0))
      {
        metric("phg_configured_t0", side, configuredT0, "PHGarfield configured t0 [bins], not a fitted edge");
      }
      if (!d.time[side])
      {
        continue;
      }
      const std::string sideName = side == 0 ? "south" : "north";
      folder(dir, "drift")->cd();
      const std::string weight = d.source[side] == "production_thresholded_sample_counts"
                                     ? "thresholded sample counts"
                                     : "sum cleaned ADC";
      if (d.rawTime[side])
      {
        d.rawTime[side]->SetOption("HIST");
        d.rawTime[side]->SetTitle(
            (plain(sideName) + ": summed raw waveform time;raw waveform sample index;" + weight).c_str());
        d.rawTime[side]->Write(("time_raw_" + sideName).c_str());
      }
      d.time[side]->SetTitle(
          (plain(sideName) + ": summed unpacked-coordinate time;unpacked TPC time bin;" + weight).c_str());
      d.time[side]->SetOption("HIST");
      d.time[side]->Write(("time_unpacked_" + sideName).c_str());
      folder(dir, "metadata")->cd();
      TNamed source(("drift_source_side" + std::to_string(side)).c_str(), d.source[side].c_str());
      source.Write();
      std::ostringstream conventionText;
      if (d.source[side] == "production_thresholded_sample_counts")
      {
        conventionText << "unpacked_time = raw_sample - " << d.appliedShift[side]
                       << "; presample shift + unpacker t0; negative transformed bins retained; raw and "
                          "unpacked spectra both saved";
      }
      else
      {
        conventionText << "Input is already unpacked TRKR_HITSET time; no presample correction applied";
      }
      TNamed shiftInfo(("time_coordinate_" + sideName).c_str(), conventionText.str().c_str());
      shiftInfo.Write();
      TH1D shiftValue(("raw_to_unpacked_shift_" + sideName).c_str(), "Applied raw-to-unpacked time shift;;bins",
                      1, .5, 1.5);
      shiftValue.SetDirectory(nullptr);
      shiftValue.SetBinContent(1, d.appliedShift[side]);
      shiftValue.Write();
      const auto data = xy(d.time[side].get(), side, false);
      const Edge t0 = edge(data.first, data.second, 2.5, 5.5, 12., 16., 5., 14., 8., 1);
      const Edge cm = edge(data.first, data.second, 180., 225., 300., 360., 210., 290., 250., -1);
      auto cleaned = data.second;
      for (size_t i = 1; i + 1 < cleaned.size(); ++i)
      {
        cleaned[i] = median({data.second[i - 1], data.second[i], data.second[i + 1]});
      }
      profile("t0", side, data.first, cleaned, t0, 2., 16., "unpacked TPC time bin");
      profile("cm_time", side, data.first, cleaned, cm, 210., 290., "unpacked TPC time bin");
      const double dt = (cm.x50 - t0.x50) * d.clock * .001;
      const bool valid = t0.valid && cm.valid && dt > 0.;
      double vd = std::numeric_limits<double>::quiet_NaN();
      double vp = vd;
      double zpad = vd;
      double zcm = vd;
      if (t0.valid)
      {
        metric("t0_bins", side, t0.x50, "t_{0} [unpacked time bins]");
      }
      if (cm.valid)
      {
        metric("cm_bins", side, cm.x50, "t_{CM} [unpacked time bins]");
      }
      if (valid)
      {
        vd = (d.pad - d.cmHalf) / dt;
        metric("data_velocity", side, vd, "v_{data} [cm/#mus]", "ratio");
      }
      if (d.z[side])
      {
        folder(dir, "drift")->cd();
        d.z[side]->SetTitle(
            (plain(sideName) + ": summed PHGarfield z;PHGarfield z [cm];sum cleaned ADC").c_str());
        d.z[side]->Write(("phgarfield_z_" + sideName).c_str());
        const auto upper = physicalZ(d.z[side].get(), side, d.pad);
        const Edge cz = edge(upper.first, upper.second, -15., -5., 5., 15., -8., 8., 0., 1, false);
        auto sub = upper.second;
        // Notebook pad background: slope from u<0, level anchored at 103..105 cm.
        TGraphErrors negative;
        for (size_t i = 0; i < upper.first.size(); ++i)
        {
          if (upper.first[i] >= -25. && upper.first[i] <= -2.)
          {
            int j = negative.GetN();
            negative.SetPoint(j, upper.first[i], upper.second[i]);
          }
        }
        bool bgValid = false;
        if (negative.GetN() >= 4)
        {
          TF1 background("dv_background", "pol1", -25., -2.);
          auto fit = negative.Fit(&background, "QSN0R");
          if (int(fit) == 0)
          {
            const double slope = background.GetParameter(1);
            std::vector<double> anchors;
            for (size_t i = 0; i < upper.first.size(); ++i)
            {
              if (upper.first[i] >= 103. && upper.first[i] <= 105.)
              {
                anchors.push_back(upper.second[i] - slope * upper.first[i]);
              }
            }
            const double intercept = median(anchors);
            if (std::isfinite(intercept))
            {
              bgValid = true;
              for (size_t i = 0; i < sub.size(); ++i)
              {
                sub[i] -= slope * upper.first[i] + intercept;
              }
            }
          }
        }
        const Edge pz =
            bgValid ? edge(upper.first, sub, 97., 100.5, 103., 105., 100., 104.9, d.pad, -1, false) : Edge{};
        profile("cm_z", side, upper.first, upper.second, cz, -15., 15., "u=sign(side)*z [cm]");
        profile("pad_z", side, upper.first, sub, pz, 97., 108., "u=sign(side)*z [cm]");
        Curve rawZ = makeCurve(sample.label + ": folded PHGarfield z upper branch", "u=sign(side)*z [cm]",
                               "sum cleaned ADC");
        Curve subZ = rawZ;
        subZ.yTitle = "background-subtracted sum ADC";
        for (size_t i = 0; i < upper.first.size(); ++i)
        {
          point(rawZ, upper.first[i], 0., upper.second[i], 0.);
          if (bgValid)
          {
            point(subZ, upper.first[i], 0., sub[i], 0.);
          }
        }
        folder(dir, "diagnostics/drift")->cd();
        graph(rawZ, "phgarfield_z_upper_branch_" + sideName)->Write();
        if (bgValid)
        {
          graph(subZ, "phgarfield_z_background_subtracted_" + sideName)->Write();
        }
        if (cz.valid)
        {
          zcm = cz.x50;
          metric("phg_cm_z", side, zcm, "PHGarfield CM edge u [cm]");
        }
        if (pz.valid)
        {
          zpad = pz.x50;
          metric("phg_pad_z", side, zpad, "PHGarfield pad edge u [cm]");
        }
        if (valid && cz.valid && pz.valid && zpad > zcm)
        {
          vp = (zpad - zcm) / dt;
          metric("phg_edge_velocity", side, vp, "effective PHGarfield mapping v [cm/#mus]", "ratio");
          metric("phg_over_data", side, vp / vd, "v_{PHG,edges}/v_{data}", "ratio");
        }
      }
      const double baseline = d.calibrations[side] > 0 ? d.phgMean[side] / d.calibrations[side]
                                                       : std::numeric_limits<double>::quiet_NaN();
      const double weighted = d.weightedADC[side] > 0 ? d.weightedV[side] / d.weightedADC[side]
                                                      : std::numeric_limits<double>::quiet_NaN();
      if (std::isfinite(baseline))
      {
        if (valid)
        {
          metric("data_over_baseline", side, vd / baseline, "v_{data}/v_{PHG,baseline}", "ratio");
        }
      }
      if (std::isfinite(weighted))
      {
        metric("phg_adc_weighted_velocity", side, weighted, "ADC-weighted PHGarfield mapping v [cm/#mus]",
               "ratio");
      }
      if (valid && std::isfinite(d.referenceV[side]))
      {
        metric("data_over_reference", side, vd / d.referenceV[side],
               "v_{data}/v_{PHG,reference}", "ratio");
      }
      csv << sample.label << ',' << side << ',' << d.source[side] << ',' << t0.valid << ',' << cm.valid << ','
          << t0.x50 << ',' << cm.x50 << ',' << std::abs(t0.x90 - t0.x10) << ',' << std::abs(cm.x90 - cm.x10)
          << ',' << vd << ',' << vp << ',' << baseline << ',' << weighted << ',' << zpad << ',' << zcm << ','
          << d.pad << ',' << d.cmHalf << ',' << d.clock << ',' << d.appliedShift[side] << '\n';
    }
    for (const auto &entry : curves)
    {
      writeCurve(entry.second, entry.first, dir);
    }
    return curves;
  }
}  // namespace TpcCalibrationQADetail
