#include "Tpc_ModuleTrackDisplay.h"

#include "tpctrackreco/Tpc_FittingTools.h"
#include "tpctrackreco/IdealPadMap.h"
#include "tpctrackreco/Tpc_ModuleTrack.h"
#include "tpctrackreco/Tpc_ModuleTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH3D.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{
  int track_color(const unsigned int itrk)
  {
    static const int colors[] = {
      kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
      kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7, kTeal + 3
    };
    return colors[itrk % (sizeof(colors) / sizeof(colors[0]))];
  }

  int layer_to_module(const unsigned int layer)
  {
    if (layer < 7 || layer > 54) return -1;
    const int module = static_cast<int>((layer - 7) / 16);
    return (module >= 0 && module < 3) ? module : -1;
  }

  unsigned long long make_unique_hit_id(const TrkrDefs::hitsetkey hsk,
                                        const TrkrDefs::hitkey hk)
  {
    return (static_cast<unsigned long long>(hsk) << 32) |
           static_cast<unsigned long long>(hk);
  }

  bool is_good_number(const double x)
  {
    return std::isfinite(x) && std::fabs(x) < 1.0e30;
  }


  double unwrap_phi_near(const double phi, const double reference)
  {
    double out = phi;
    while (out - reference > TMath::Pi()) out -= 2.0 * TMath::Pi();
    while (out - reference < -TMath::Pi()) out += 2.0 * TMath::Pi();
    return out;
  }

  bool sagitta_y_at_x(const double x,
                      const double mean_y,
                      const Tpc_FittingTools::SagittaFit& fit,
                      double& y_out)
  {
    if (!fit.ok || !is_good_number(x) || !is_good_number(mean_y)) return false;
    if (!is_good_number(fit.S) || !is_good_number(fit.x0) ||
        !is_good_number(fit.invR) || !is_good_number(fit.theta) ||
        !is_good_number(fit.b)) return false;

    const double c = std::cos(fit.theta);
    const double st = std::sin(fit.theta);
    const double ymin = mean_y - 0.8;
    const double ymax = mean_y + 0.8;
    const int nscan = 240;

    bool have_best = false;
    double best_y = 0.0;
    double best_dist = 1.0e99;

    double yprev = ymin;
    double xprev_rot = c * x + st * (yprev - fit.b);
    double yprev_rot = -st * x + c * (yprev - fit.b);
    double gprev = yprev_rot - Tpc_FittingTools::sagittaModel(xprev_rot, fit.S, fit.x0, fit.invR);

    double best_abs_g = std::fabs(gprev);
    double best_abs_y = yprev;

    for (int iscan = 1; iscan <= nscan; ++iscan)
    {
      const double ycur = ymin + (ymax - ymin) * static_cast<double>(iscan) / static_cast<double>(nscan);
      const double xcur_rot = c * x + st * (ycur - fit.b);
      const double ycur_rot = -st * x + c * (ycur - fit.b);
      const double gcur = ycur_rot - Tpc_FittingTools::sagittaModel(xcur_rot, fit.S, fit.x0, fit.invR);

      if (std::fabs(gcur) < best_abs_g)
      {
        best_abs_g = std::fabs(gcur);
        best_abs_y = ycur;
      }

      const bool bracket =
        gprev == 0.0 || gcur == 0.0 ||
        (gprev < 0.0 && gcur > 0.0) ||
        (gprev > 0.0 && gcur < 0.0);

      if (bracket)
      {
        double ya = yprev;
        double yb = ycur;
        double ga = gprev;

        for (int ib = 0; ib < 50; ++ib)
        {
          const double ym = 0.5 * (ya + yb);
          const double xm_rot = c * x + st * (ym - fit.b);
          const double ym_rot = -st * x + c * (ym - fit.b);
          const double gm = ym_rot - Tpc_FittingTools::sagittaModel(xm_rot, fit.S, fit.x0, fit.invR);

          if ((ga < 0.0 && gm <= 0.0) || (ga > 0.0 && gm >= 0.0))
          {
            ya = ym;
            ga = gm;
          }
          else
          {
            yb = ym;
          }
        }

        const double y = 0.5 * (ya + yb);
        const double dist = std::fabs(y - mean_y);
        if (!have_best || dist < best_dist)
        {
          have_best = true;
          best_dist = dist;
          best_y = y;
        }
      }

      yprev = ycur;
      gprev = gcur;
    }

    if (have_best)
    {
      y_out = best_y;
      return true;
    }

    if (best_abs_g < 2.0e-2)
    {
      y_out = best_abs_y;
      return true;
    }

    return false;
  }

  void style_hit_graph(TGraph* g, const int color)
  {
    if (!g) return;
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
  }

  void style_fit_graph(TGraph* g, const int color)
  {
    if (!g) return;
    g->SetMarkerStyle(1);
    g->SetMarkerSize(0.0);
    g->SetLineColor(color);
    g->SetLineWidth(3);
    g->SetLineStyle(2);
  }

  void style_fit_line_3d(TPolyLine3D* line, const int color)
  {
    if (!line) return;
    line->SetLineColor(color);
    line->SetLineWidth(4);
    line->SetLineStyle(1);
  }

  struct GroupKey
  {
    GroupKey() : side(-1), sector(-1), module(-1) {}
    GroupKey(const int s, const int sec, const int mod)
      : side(s), sector(sec), module(mod) {}

    int side;
    int sector;
    int module;

    bool operator<(const GroupKey& rhs) const
    {
      if (side != rhs.side) return side < rhs.side;
      if (sector != rhs.sector) return sector < rhs.sector;
      return module < rhs.module;
    }
  };

  struct FitResult
  {
    FitResult()
      : ok(false)
      , use_sagitta(false)
      , phi_slope(0.0)
      , phi_intercept(0.0)
      , tbin_slope(0.0)
      , tbin_intercept(0.0)
      , phi_sagitta()
      , mean_phi(0.0)
      , rmin(0.0)
      , rmax(0.0)
    {}

    bool ok;
    bool use_sagitta;
    double phi_slope;
    double phi_intercept;
    double tbin_slope;
    double tbin_intercept;
    Tpc_FittingTools::SagittaFit phi_sagitta;
    double mean_phi;
    double rmin;
    double rmax;
  };

  struct HardwareFitResult
  {
    HardwareFitResult()
      : ok(false)
      , use_sagitta(false)
      , pad_slope(0.0)
      , pad_intercept(0.0)
      , tbin_slope(0.0)
      , tbin_intercept(0.0)
      , pad_sagitta()
      , mean_pad(0.0)
      , layer_min(0.0)
      , layer_max(0.0)
    {}

    bool ok;
    bool use_sagitta;
    double pad_slope;
    double pad_intercept;
    double tbin_slope;
    double tbin_intercept;
    Tpc_FittingTools::SagittaFit pad_sagitta;
    double mean_pad;
    double layer_min;
    double layer_max;
  };

  struct GraphBundle
  {
    GraphBundle()
      : mg_phi_radius_hits(nullptr)
      , mg_tbin_radius_hits(nullptr)
      , mg_tbin_phi_hits(nullptr)
      , mg_phi_radius_fits(nullptr)
      , mg_tbin_radius_fits(nullptr)
      , mg_tbin_phi_fits(nullptr)
      , h3_hardware_hits(nullptr)
      , h3_adc_hits(nullptr)
      , h3_adc_unassociated_hits(nullptr)
      , c3_hardware_hits_fits(nullptr)
      , c3_adc_hits_fits(nullptr)
      , c3_adc_unassociated_hits(nullptr)
      , phi_reference(0.0)
    {}

    TMultiGraph* mg_phi_radius_hits;
    TMultiGraph* mg_tbin_radius_hits;
    TMultiGraph* mg_tbin_phi_hits;

    TMultiGraph* mg_phi_radius_fits;
    TMultiGraph* mg_tbin_radius_fits;
    TMultiGraph* mg_tbin_phi_fits;

    // x = timebin, y = pad, z = layer, weight = ADC
    TH3D* h3_hardware_hits;

    // x = timebin, y = phi, z = radius, weight = ADC
    TH3D* h3_adc_hits;
    TH3D* h3_adc_unassociated_hits;

    TCanvas* c3_hardware_hits_fits;
    TCanvas* c3_adc_hits_fits;
    TCanvas* c3_adc_unassociated_hits;

    double phi_reference;

    std::vector<TPolyLine3D*> hardware_fit_lines_3d;
    std::vector<std::string> hardware_fit_line_names_3d;

    std::vector<TPolyLine3D*> fit_lines_3d;
    std::vector<std::string> fit_line_names_3d;

    std::set<unsigned long long> filled_hit_ids;
  };

  void make_bundle_graphs(GraphBundle& b,
                          const GroupKey& key,
                          const unsigned int evt,
                          const IdealPadMap* idealPadMap)
  {
    if (!idealPadMap) return;

    const TString tag = Form("s%d_sec%02d_mod%d", key.side, key.sector, key.module);

    b.mg_phi_radius_hits = new TMultiGraph();
    b.mg_phi_radius_hits->SetName(Form("mg_%s_phi_radius_hits", tag.Data()));
    b.mg_phi_radius_hits->SetTitle(Form("event %u side %d sector %d module %d hits;#phi;radius [cm]",
                                        evt, key.side, key.sector, key.module));

    b.mg_tbin_radius_hits = new TMultiGraph();
    b.mg_tbin_radius_hits->SetName(Form("mg_%s_tbin_radius_hits", tag.Data()));
    b.mg_tbin_radius_hits->SetTitle(Form("event %u side %d sector %d module %d hits;timebin;radius [cm]",
                                         evt, key.side, key.sector, key.module));

    b.mg_tbin_phi_hits = new TMultiGraph();
    b.mg_tbin_phi_hits->SetName(Form("mg_%s_tbin_phi_hits", tag.Data()));
    b.mg_tbin_phi_hits->SetTitle(Form("event %u side %d sector %d module %d hits;timebin;#phi",
                                      evt, key.side, key.sector, key.module));

    b.mg_phi_radius_fits = new TMultiGraph();
    b.mg_phi_radius_fits->SetName(Form("mg_%s_phi_radius_fits", tag.Data()));
    b.mg_phi_radius_fits->SetTitle(Form("event %u side %d sector %d module %d display fits;#phi;radius [cm]",
                                        evt, key.side, key.sector, key.module));

    b.mg_tbin_radius_fits = new TMultiGraph();
    b.mg_tbin_radius_fits->SetName(Form("mg_%s_tbin_radius_fits", tag.Data()));
    b.mg_tbin_radius_fits->SetTitle(Form("event %u side %d sector %d module %d display fits;timebin;radius [cm]",
                                         evt, key.side, key.sector, key.module));

    b.mg_tbin_phi_fits = new TMultiGraph();
    b.mg_tbin_phi_fits->SetName(Form("mg_%s_tbin_phi_fits", tag.Data()));
    b.mg_tbin_phi_fits->SetTitle(Form("event %u side %d sector %d module %d display fits;timebin;#phi",
                                      evt, key.side, key.sector, key.module));

    const int layer_min = 7 + 16 * key.module;
    const int layer_max = layer_min + 15;
    const double radius_min = idealPadMap->get_radius(static_cast<unsigned int>(layer_min));
    const double radius_max = idealPadMap->get_radius(static_cast<unsigned int>(layer_max));
    const double radius_min_guess = std::min(radius_min, radius_max) - 1.0;
    const double radius_max_guess = std::max(radius_min, radius_max) + 1.0;

    const unsigned int pads_per_sector =
      idealPadMap->get_pads_per_sector(static_cast<unsigned int>(key.module));

    if (pads_per_sector < 2U) return;

    const unsigned int nPads = pads_per_sector;
    const int pad_min = static_cast<int>(nPads) * key.sector;
    const int pad_max = pad_min + static_cast<int>(nPads) - 1;
    const unsigned int layer_ref = static_cast<unsigned int>(layer_min);

    // first and last pad centers in this sector
    const double phi_first =
      idealPadMap->get_phi(key.side,
                           static_cast<unsigned int>(key.sector),
                           layer_ref,
                           0);

    const double phi_last_wrapped =
      idealPadMap->get_phi(key.side,
                           static_cast<unsigned int>(key.sector),
                           layer_ref,
                           pads_per_sector - 1);

    const double phi_last = unwrap_phi_near(phi_last_wrapped, phi_first);
    const double dphi = (phi_last - phi_first) / static_cast<double>(pads_per_sector - 1);

    double phi_min = phi_first - 0.5 * dphi;
    double phi_max = phi_last + 0.5 * dphi;
    if (phi_max < phi_min) std::swap(phi_min, phi_max);
    b.phi_reference = 0.5 * (phi_min + phi_max);

    b.h3_hardware_hits = new TH3D(Form("h3_%s_hardware_hits", tag.Data()),
                                  Form("event %u side %d sector %d module %d associated hits;timebin;pad;layer",
                                       evt, key.side, key.sector, key.module),
                                  512, -0.5, 511.5,
                                  static_cast<int>(nPads), static_cast<double>(pad_min) - 0.5, static_cast<double>(pad_max) + 0.5,
                                  16, static_cast<double>(layer_min) - 0.5, static_cast<double>(layer_max) + 0.5);
    b.h3_hardware_hits->SetStats(0);

    b.h3_adc_hits = new TH3D(Form("h3_%s_adc_hits", tag.Data()),
                             Form("event %u side %d sector %d module %d associated hits;timebin;#phi;radius [cm]",
                                  evt, key.side, key.sector, key.module),
                             512, -0.5, 511.5,
                             static_cast<int>(nPads), phi_min, phi_max,
                             20, radius_min_guess, radius_max_guess);
    b.h3_adc_hits->SetStats(0);

    b.h3_adc_unassociated_hits = new TH3D(Form("h3_%s_adc_unassociated_hits", tag.Data()),
              Form("event %u side %d sector %d module %d unassociated hits;timebin;#phi;radius [cm]",
                    evt, key.side, key.sector, key.module),
              512, -0.5, 511.5,
              static_cast<int>(nPads), phi_min, phi_max,
              16, radius_min_guess, radius_max_guess);
    b.h3_adc_unassociated_hits->SetStats(0);

  }

  GraphBundle& get_bundle(std::map<GroupKey, GraphBundle>& bundles,
                          const GroupKey& key,
                          const unsigned int evt,
                          const IdealPadMap* idealPadMap)
  {
    std::map<GroupKey, GraphBundle>::iterator it = bundles.find(key);
    if (it == bundles.end())
    {
      const std::pair<std::map<GroupKey, GraphBundle>::iterator, bool> inserted =
        bundles.insert(std::make_pair(key, GraphBundle()));
      it = inserted.first;
      make_bundle_graphs(it->second, key, evt, idealPadMap);
    }
    return it->second;
  }

  HardwareFitResult fit_hardware_track_points(const std::vector<Tpc_ModuleTrackDisplay::HitPoint>& pts,
                                               const int fit_mode,
                                               const double weight_power,
                                               const double weight_floor_frac)
  {
    HardwareFitResult result;
    if (pts.size() < 2) return result;

    double max_adc = 0.0;
    for (const auto& p : pts) max_adc = std::max(max_adc, static_cast<double>(p.adc));
    if (max_adc <= 0.0) max_adc = 1.0;

    std::vector<Tpc_FittingTools::FitPoint> layer_pad_points;
    std::vector<Tpc_FittingTools::FitPoint> layer_tbin_points;
    layer_pad_points.reserve(pts.size());
    layer_tbin_points.reserve(pts.size());

    result.layer_min = static_cast<double>(pts.front().layer);
    result.layer_max = static_cast<double>(pts.front().layer);
    result.mean_pad = 0.0;

    for (const auto& p : pts)
    {
      const double layer = static_cast<double>(p.layer);
      const double pad = static_cast<double>(p.pad);
      const double tbin = static_cast<double>(p.tbin);
      const double w = Tpc_FittingTools::adcWeight(static_cast<double>(p.adc), max_adc,
                                         weight_power, weight_floor_frac);

      layer_pad_points.emplace_back(layer, pad, w);
      layer_tbin_points.emplace_back(layer, tbin, w);
      result.mean_pad += pad;
      result.layer_min = std::min(result.layer_min, layer);
      result.layer_max = std::max(result.layer_max, layer);
    }

    result.mean_pad /= static_cast<double>(layer_pad_points.size());

    const Tpc_FittingTools::LineFit pad_fit = Tpc_FittingTools::fitLine(layer_pad_points);
    const Tpc_FittingTools::LineFit tbin_fit = Tpc_FittingTools::fitLine(layer_tbin_points);

    if (!pad_fit.ok || !tbin_fit.ok) return result;

    result.ok = true;
    result.pad_slope = pad_fit.slope;
    result.pad_intercept = pad_fit.intercept;
    result.tbin_slope = tbin_fit.slope;
    result.tbin_intercept = tbin_fit.intercept;

    if (fit_mode == Tpc_FittingTools::FIT_SAGITTA)
    {
      result.pad_sagitta = Tpc_FittingTools::fitSagitta(layer_pad_points);
      result.use_sagitta = result.pad_sagitta.ok;
    }

    return result;
  }

  FitResult fit_track_points(const std::vector<Tpc_ModuleTrackDisplay::HitPoint>& pts,
                             const int fit_mode,
                             const double weight_power,
                             const double weight_floor_frac)
  {
    FitResult result;
    if (pts.size() < 2) return result;

    double max_adc = 0.0;
    for (const auto& p : pts) max_adc = std::max(max_adc, static_cast<double>(p.adc));
    if (max_adc <= 0.0) max_adc = 1.0;

    const double phi_ref = pts.front().phi;
    std::vector<Tpc_FittingTools::FitPoint> radius_phi_points;
    std::vector<Tpc_FittingTools::FitPoint> radius_tbin_points;
    radius_phi_points.reserve(pts.size());
    radius_tbin_points.reserve(pts.size());

    result.rmin = pts.front().radius;
    result.rmax = pts.front().radius;
    result.mean_phi = 0.0;

    for (const auto& p : pts)
    {
      if (!is_good_number(p.radius) || !is_good_number(p.phi)) continue;
      const double w = Tpc_FittingTools::adcWeight(static_cast<double>(p.adc), max_adc,
                                         weight_power, weight_floor_frac);
      const double phi_unwrapped = unwrap_phi_near(p.phi, phi_ref);
      radius_phi_points.emplace_back(p.radius, phi_unwrapped, w);
      radius_tbin_points.emplace_back(p.radius, static_cast<double>(p.tbin), w);
      result.mean_phi += phi_unwrapped;
      result.rmin = std::min(result.rmin, p.radius);
      result.rmax = std::max(result.rmax, p.radius);
    }

    if (radius_phi_points.empty()) return result;
    result.mean_phi /= static_cast<double>(radius_phi_points.size());

    const Tpc_FittingTools::LineFit phi_fit = Tpc_FittingTools::fitLine(radius_phi_points);
    const Tpc_FittingTools::LineFit tbin_fit = Tpc_FittingTools::fitLine(radius_tbin_points);

    if (!phi_fit.ok || !tbin_fit.ok) return result;

    result.ok = true;
    result.phi_slope = phi_fit.slope;
    result.phi_intercept = phi_fit.intercept;
    result.tbin_slope = tbin_fit.slope;
    result.tbin_intercept = tbin_fit.intercept;

    if (fit_mode == Tpc_FittingTools::FIT_SAGITTA)
    {
      result.phi_sagitta = Tpc_FittingTools::fitSagitta(radius_phi_points);
      result.use_sagitta = result.phi_sagitta.ok;
    }

    return result;
  }

  void add_fit_objects(GraphBundle& b,
                       const GroupKey& key,
                       const unsigned int evt,
                       const unsigned int tid,
                       const int color,
                       const FitResult& fit)
  {
    if (!fit.ok || fit.rmax <= fit.rmin) return;

    const TString tag = Form("s%d_sec%02d_mod%d_trk%u", key.side, key.sector, key.module, tid);

    TGraph* g_phi_radius_fit = new TGraph();
    g_phi_radius_fit->SetName(Form("g_%s_phi_radius_fit", tag.Data()));
    g_phi_radius_fit->SetTitle(Form("event %u track %u fit;#phi;radius [cm]", evt, tid));
    style_fit_graph(g_phi_radius_fit, color);

    TGraph* g_tbin_radius_fit = new TGraph();
    g_tbin_radius_fit->SetName(Form("g_%s_tbin_radius_fit", tag.Data()));
    g_tbin_radius_fit->SetTitle(Form("event %u track %u fit;timebin;radius [cm]", evt, tid));
    style_fit_graph(g_tbin_radius_fit, color);

    TGraph* g_tbin_phi_fit = new TGraph();
    g_tbin_phi_fit->SetName(Form("g_%s_tbin_phi_fit", tag.Data()));
    g_tbin_phi_fit->SetTitle(Form("event %u track %u fit;timebin;#phi", evt, tid));
    style_fit_graph(g_tbin_phi_fit, color);

    const int npts = 51;
    TPolyLine3D* line3 = new TPolyLine3D(npts);
    const std::string line_name = Form("line3_%s_tbin_phi_radius_fit", tag.Data());

    for (int i = 0; i < npts; ++i)
    {
      const double f = static_cast<double>(i) / static_cast<double>(npts - 1);
      const double r = fit.rmin + f * (fit.rmax - fit.rmin);
      double phi = fit.phi_slope * r + fit.phi_intercept;
      if (fit.use_sagitta)
      {
        double phi_sagitta = 0.0;
        if (sagitta_y_at_x(r, fit.mean_phi, fit.phi_sagitta, phi_sagitta)) phi = phi_sagitta;
      }
      const double tbin = fit.tbin_slope * r + fit.tbin_intercept;
      const double phi_draw = unwrap_phi_near(phi, b.phi_reference);

      g_phi_radius_fit->SetPoint(i, phi_draw, r);
      g_tbin_radius_fit->SetPoint(i, tbin, r);
      g_tbin_phi_fit->SetPoint(i, tbin, phi_draw);
      line3->SetPoint(i, tbin, phi_draw, r);
    }

    style_fit_line_3d(line3, color);
    b.fit_lines_3d.push_back(line3);
    b.fit_line_names_3d.push_back(line_name);

    b.mg_phi_radius_fits->Add(g_phi_radius_fit, "L");
    b.mg_tbin_radius_fits->Add(g_tbin_radius_fit, "L");
    b.mg_tbin_phi_fits->Add(g_tbin_phi_fit, "L");
  }

  void add_hardware_fit_objects(GraphBundle& b,
                                const GroupKey& key,
                                const unsigned int tid,
                                const int color,
                                const HardwareFitResult& fit)
  {
    if (!fit.ok || fit.layer_max <= fit.layer_min) return;

    const TString tag = Form("s%d_sec%02d_mod%d_trk%u", key.side, key.sector, key.module, tid);

    const int npts = 51;
    TPolyLine3D* line3 = new TPolyLine3D(npts);
    const std::string line_name = Form("line3_%s_layer_pad_tbin_fit", tag.Data());

    for (int i = 0; i < npts; ++i)
    {
      const double f = static_cast<double>(i) / static_cast<double>(npts - 1);
      const double layer = fit.layer_min + f * (fit.layer_max - fit.layer_min);
      double pad = fit.pad_slope * layer + fit.pad_intercept;
      if (fit.use_sagitta)
      {
        double pad_sagitta = 0.0;
        if (sagitta_y_at_x(layer, fit.mean_pad, fit.pad_sagitta, pad_sagitta)) pad = pad_sagitta;
      }
      const double tbin = fit.tbin_slope * layer + fit.tbin_intercept;
      line3->SetPoint(i, tbin, pad, layer);
    }

    style_fit_line_3d(line3, color);
    b.hardware_fit_lines_3d.push_back(line3);
    b.hardware_fit_line_names_3d.push_back(line_name);
  }

  void write_bundle(GraphBundle& b, const GroupKey& key)
  {
    const TString tag = Form("s%d_sec%02d_mod%d", key.side, key.sector, key.module);

    if (b.h3_hardware_hits)
    {
      b.c3_hardware_hits_fits = new TCanvas(Form("c3_%s_hardware_hits_fits", tag.Data()),
                                            Form("side %d sector %d module %d hardware hits and display fits",
                                                 key.side, key.sector, key.module),
                                            1200, 900);
      b.h3_hardware_hits->Draw("BOX2Z");
      for (unsigned int i = 0; i < b.hardware_fit_lines_3d.size(); ++i)
      {
        if (!b.hardware_fit_lines_3d[i]) continue;
        b.hardware_fit_lines_3d[i]->Draw("same");
      }
      b.c3_hardware_hits_fits->Modified();
      b.c3_hardware_hits_fits->Update();
      b.c3_hardware_hits_fits->Write();
    }

    if (b.h3_adc_hits)
    {
      b.c3_adc_hits_fits = new TCanvas(Form("c3_%s_adc_hits_fits", tag.Data()),
                                       Form("side %d sector %d module %d ADC hits and display fits",
                                            key.side, key.sector, key.module),
                                       1200, 900);
      b.h3_adc_hits->Draw("BOX2Z");
      for (unsigned int i = 0; i < b.fit_lines_3d.size(); ++i)
      {
        if (!b.fit_lines_3d[i]) continue;
        b.fit_lines_3d[i]->Draw("same");
        //b.fit_lines_3d[i]->Write(b.fit_line_names_3d[i].c_str());
      }
      b.c3_adc_hits_fits->Modified();
      b.c3_adc_hits_fits->Update();
      b.c3_adc_hits_fits->Write();
    }

    if (b.h3_adc_unassociated_hits)
    {
      b.c3_adc_unassociated_hits = new TCanvas(Form("c3_%s_adc_unassociated_hits", tag.Data()),
                                               Form("side %d sector %d module %d unassociated ADC hits",
                                                    key.side, key.sector, key.module),
                                               1200, 900);
      b.h3_adc_unassociated_hits->Draw("BOX2Z");
      b.c3_adc_unassociated_hits->Modified();
      b.c3_adc_unassociated_hits->Update();
      b.c3_adc_unassociated_hits->Write();
    }
/*
    if (b.mg_phi_radius_hits) b.mg_phi_radius_hits->Write();
    if (b.mg_tbin_radius_hits) b.mg_tbin_radius_hits->Write();
    if (b.mg_tbin_phi_hits) b.mg_tbin_phi_hits->Write();
    if (b.mg_phi_radius_fits) b.mg_phi_radius_fits->Write();
    if (b.mg_tbin_radius_fits) b.mg_tbin_radius_fits->Write();
    if (b.mg_tbin_phi_fits) b.mg_tbin_phi_fits->Write();
    if (b.h3_adc_hits) b.h3_adc_hits->Write();
    if (b.h3_adc_unassociated_hits) b.h3_adc_unassociated_hits->Write();
    */
  }
}

Tpc_ModuleTrackDisplay::HitPoint::HitPoint()
  : ok(false)
  , hitsetkey(0)
  , hitkey(0)
  , side(-1)
  , sector(0)
  , layer(0)
  , pad(0)
  , tbin(0)
  , adc(0)
  , radius(0.0)
  , phi(0.0)
{}

Tpc_ModuleTrackDisplay::Tpc_ModuleTrackDisplay(const std::string& name,
                                           const std::string& outfilename,
                                           const std::string& trackNodeName,
                                           unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_trackNodeName(trackNodeName)
  , m_maxEventDisplays(maxEventDisplays)
  , m_evt(0)
  , m_eventsSaved(0)
  , m_outfile(nullptr)
  , m_tracks(nullptr)
  , m_hits(nullptr)
  , m_idealPadMap(new IdealPadMap())
  , m_fitMode(Tpc_FittingTools::FIT_SAGITTA)
  , m_fitWeightPower(1.0)
  , m_fitWeightFloorFrac(0.05)
{}

Tpc_ModuleTrackDisplay::~Tpc_ModuleTrackDisplay()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
}

int Tpc_ModuleTrackDisplay::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << "Tpc_ModuleTrackDisplay::Init - cannot open output file "
              << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_idealPadMap)
  {
    m_idealPadMap = new IdealPadMap();
  }

  if (!m_idealPadMap->is_loaded() && m_idealPadMap->load_from_cdb(Verbosity()) != 0)
  {
    std::cerr << "Tpc_ModuleTrackDisplay::Init - failed to load IdealPadMap from CDB"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_outfile->mkdir("events");
  std::cout << "Tpc_ModuleTrackDisplay::Init - writing display fits to "
            << m_outfilename << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_ModuleTrackDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_evt;

  if (!get_nodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;

  const unsigned int ntracks = m_tracks ? m_tracks->size() : 0;
  if (m_eventsSaved >= m_maxEventDisplays) return Fun4AllReturnCodes::EVENT_OK;

  TDirectory* eventsTop = m_outfile->GetDirectory("events");
  if (!eventsTop) eventsTop = m_outfile->mkdir("events");

  eventsTop->cd();
  TDirectory* eventDir = eventsTop->mkdir(Form("event_%06u", m_evt));
  if (!eventDir)
  {
    std::cerr << "Tpc_ModuleTrackDisplay::process_event - failed to create event directory"
              << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  eventDir->cd();

  std::map<GroupKey, GraphBundle> bundles;
  std::set<unsigned long long> associated_hit_ids;

  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const Tpc_ModuleTrack* trk = m_tracks->get_track(itrk);
    if (!trk) continue;

    const unsigned int tid = trk->get_track_id();
    const int color = track_color(itrk);

    std::map<GroupKey, std::vector<HitPoint> > points_by_group;

    for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
    {
      const Tpc_ModuleTrack::HitIndex idx = trk->get_hit_index(ih);
      const HitPoint p = make_hit_point(idx.first, idx.second);
      if (!p.ok) continue;

      associated_hit_ids.insert(make_unique_hit_id(p.hitsetkey, p.hitkey));

      const int module = layer_to_module(p.layer);
      if (module < 0) continue;

      const GroupKey key(p.side, static_cast<int>(p.sector), module);
      points_by_group[key].push_back(p);
    }

    for (std::map<GroupKey, std::vector<HitPoint> >::iterator pg = points_by_group.begin();
         pg != points_by_group.end(); ++pg)
    {
      const GroupKey& key = pg->first;
      const std::vector<HitPoint>& pts = pg->second;
      if (pts.empty()) continue;

      GraphBundle& b = get_bundle(bundles, key, m_evt, m_idealPadMap);

      TGraph* g_phi_radius_hits = new TGraph();
      g_phi_radius_hits->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_phi_radius_hits",
                                      key.side, key.sector, key.module, tid));
      g_phi_radius_hits->SetTitle(Form("event %u track %u hits;#phi;radius [cm]", m_evt, tid));
      style_hit_graph(g_phi_radius_hits, color);

      TGraph* g_tbin_radius_hits = new TGraph();
      g_tbin_radius_hits->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_tbin_radius_hits",
                                       key.side, key.sector, key.module, tid));
      g_tbin_radius_hits->SetTitle(Form("event %u track %u hits;timebin;radius [cm]", m_evt, tid));
      style_hit_graph(g_tbin_radius_hits, color);

      TGraph* g_tbin_phi_hits = new TGraph();
      g_tbin_phi_hits->SetName(Form("g_s%d_sec%02d_mod%d_trk%u_tbin_phi_hits",
                                    key.side, key.sector, key.module, tid));
      g_tbin_phi_hits->SetTitle(Form("event %u track %u hits;timebin;#phi", m_evt, tid));
      style_hit_graph(g_tbin_phi_hits, color);

      for (unsigned int ip = 0; ip < pts.size(); ++ip)
      {
        const HitPoint& p = pts[ip];
        const int n = g_phi_radius_hits->GetN();
        const double phi_draw = unwrap_phi_near(p.phi, b.phi_reference);

        g_phi_radius_hits->SetPoint(n, phi_draw, p.radius);
        g_tbin_radius_hits->SetPoint(n, static_cast<double>(p.tbin), p.radius);
        g_tbin_phi_hits->SetPoint(n, static_cast<double>(p.tbin), phi_draw);

        const unsigned long long uid = make_unique_hit_id(p.hitsetkey, p.hitkey);
        if (b.filled_hit_ids.insert(uid).second)
        {
          if (b.h3_hardware_hits)
          {
            b.h3_hardware_hits->Fill(static_cast<double>(p.tbin),
                                     static_cast<double>(p.pad),
                                     static_cast<double>(p.layer),
                                     static_cast<double>(p.adc));
          }
          if (b.h3_adc_hits)
          {
            b.h3_adc_hits->Fill(static_cast<double>(p.tbin),
                                phi_draw,
                                p.radius,
                                static_cast<double>(p.adc));
          }
        }
      }

      b.mg_phi_radius_hits->Add(g_phi_radius_hits, "LP");
      b.mg_tbin_radius_hits->Add(g_tbin_radius_hits, "LP");
      b.mg_tbin_phi_hits->Add(g_tbin_phi_hits, "LP");

      const FitResult fit = fit_track_points(pts, m_fitMode, m_fitWeightPower, m_fitWeightFloorFrac);
      add_fit_objects(b, key, m_evt, tid, color, fit);

      const HardwareFitResult hardwareFit =
        fit_hardware_track_points(pts, m_fitMode, m_fitWeightPower, m_fitWeightFloorFrac);
      add_hardware_fit_objects(b, key, tid, color, hardwareFit);
    }
  }

  // Fill all TRKR_HITSET hits that are not associated with any reconstructed in-module track.
  TrkrHitSetContainer::ConstRange hitset_range = m_hits->getHitSets();
  for (TrkrHitSetContainer::ConstIterator hsiter = hitset_range.first;
       hsiter != hitset_range.second; ++hsiter)
  {
    TrkrHitSet* hitset = hsiter->second;
    if (!hitset) continue;

    const TrkrDefs::hitsetkey hsk = hsiter->first;
    const unsigned int layer = TrkrDefs::getLayer(hsk);
    const int module = layer_to_module(layer);
    if (module < 0) continue;

    const int side = static_cast<int>(TpcDefs::getSide(hsk));
    const int sector = static_cast<int>(TpcDefs::getSectorId(hsk));
    const GroupKey key(side, sector, module);
    GraphBundle& b = get_bundle(bundles, key, m_evt, m_idealPadMap);

    TrkrHitSet::ConstRange hit_range = hitset->getHits();
    for (TrkrHitSet::ConstIterator hiter = hit_range.first;
         hiter != hit_range.second; ++hiter)
    {
      const TrkrDefs::hitkey hk = hiter->first;
      const unsigned long long uid = make_unique_hit_id(hsk, hk);
      if (associated_hit_ids.find(uid) != associated_hit_ids.end()) continue;

      const HitPoint p = make_hit_point(hsk, hk);
      if (!p.ok) continue;

      if (b.h3_adc_unassociated_hits)
      {
        b.h3_adc_unassociated_hits->Fill(static_cast<double>(p.tbin),
                                         unwrap_phi_near(p.phi, b.phi_reference),
                                         p.radius,
                                         static_cast<double>(p.adc));
      }
    }
  }

  for (std::map<GroupKey, GraphBundle>::iterator ib = bundles.begin();
       ib != bundles.end(); ++ib)
  {
    eventDir->cd();
    write_bundle(ib->second, ib->first);
  }

  std::cout << "Tpc_ModuleTrackDisplay - saved event "
            << m_evt << " with " << ntracks
            << " tracks and " << bundles.size()
            << " side/sector/module groups" << std::endl;

  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_ModuleTrackDisplay::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    std::cout << "Tpc_ModuleTrackDisplay::End - events seen: " << m_evt << std::endl;
    m_outfile->cd();
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
  }

  std::cout << "Tpc_ModuleTrackDisplay::End - wrote " << m_outfilename << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool Tpc_ModuleTrackDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_tracks = findNode::getClass<Tpc_ModuleTrackContainer>(topNode, m_trackNodeName.c_str());

  if (!m_tracks)
  {
    const char* candidate_names[] = {
      "TPC_MODULETRACKS",
      "Tpc_ModuleTrackReco",
      "Tpc_ModuleTrackContainer",
      "TPC_MODULETRACKCONTAINER",
      "TPC_MODULETRACKS_CONTAINER"
    };

    for (unsigned int i = 0;
         i < sizeof(candidate_names) / sizeof(candidate_names[0]) && !m_tracks;
         ++i)
    {
      m_tracks = findNode::getClass<Tpc_ModuleTrackContainer>(topNode, candidate_names[i]);
      if (m_tracks) m_trackNodeName = candidate_names[i];
    }
  }

  // IdealPadMap is not a Fun4All node/container here.
  // The display owns it and uses it only to translate hardware coordinates
  // (side, sector, layer, pad) -> (phi, radius).
  if (!m_idealPadMap)
  {
    m_idealPadMap = new IdealPadMap();
  }

  if (!m_tracks)
  {
    std::cerr << "Tpc_ModuleTrackDisplay - could not find Tpc_ModuleTrackContainer node" << std::endl;
    return false;
  }

  if (!m_hits)
  {
    std::cerr << "Tpc_ModuleTrackDisplay - missing TRKR_HITSET" << std::endl;
    return false;
  }

  if (!m_idealPadMap)
  {
    std::cerr << "Tpc_ModuleTrackDisplay - failed to create IdealPadMap" << std::endl;
    return false;
  }

  return true;
}

Tpc_ModuleTrackDisplay::HitPoint
Tpc_ModuleTrackDisplay::make_hit_point(const TrkrDefs::hitsetkey hsk,
                                     const TrkrDefs::hitkey hk) const
{
  HitPoint p;

  TrkrHitSet* hitset = m_hits ? m_hits->findHitSet(hsk) : nullptr;
  if (!hitset) return p;

  TrkrHit* hit = hitset->getHit(hk);
  if (!hit) return p;

  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.side = static_cast<int>(TpcDefs::getSide(hsk));
  p.sector = static_cast<unsigned int>(TpcDefs::getSectorId(hsk));
  p.layer = TrkrDefs::getLayer(hsk);
  p.pad = TpcDefs::getPad(hk);
  p.tbin = TpcDefs::getTBin(hk);
  p.adc = hit->getAdc();

  if (layer_to_module(p.layer) < 0) return p;

  p.radius = ideal_radius(p.layer);
  p.phi = ideal_phi(p.side, p.sector, p.layer, p.pad);

  if (!is_good_number(p.radius) || !is_good_number(p.phi)) return p;

  p.ok = true;
  return p;
}

double Tpc_ModuleTrackDisplay::ideal_radius(const unsigned int layer) const
{
  if (!m_idealPadMap) return 0.0;

  // IdealPadMap is intentionally the only geometry source used by this display.
  // If your IdealPadMap method names differ, this is the only radius call to adjust.
  return m_idealPadMap->get_radius(layer);
}

double Tpc_ModuleTrackDisplay::ideal_phi(const int side,
                                       const unsigned int sector,
                                       const unsigned int layer,
                                       const unsigned int pad) const
{
  if (!m_idealPadMap) return 0.0;

  // Keep HitPoint::pad in raw hardware coordinates for display, but IdealPadMap
  // expects the pad index local to the sector.
  const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(layer);
  if (pads_per_sector == 0U) return 0.0;
  const unsigned int local_pad = pad % pads_per_sector;
  return m_idealPadMap->get_phi(side, sector, layer, local_pad);
}
