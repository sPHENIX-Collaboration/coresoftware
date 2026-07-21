#include "Tpc_AssembledTrackDisplay.h"

#include "tpctrackreco/Tpc_FittingTools.h"
#include "tpctrackreco/Tpc_AssembledTrack.h"
#include "tpctrackreco/Tpc_AssembledTrackContainer.h"
#include "tpctrackreco/IdealPadMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TMath.h>
#include <TPolyLine3D.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    while (phi > TMath::Pi()) phi -= 2.0 * TMath::Pi();
    while (phi <= -TMath::Pi()) phi += 2.0 * TMath::Pi();
    return phi;
  }

  double unwrap_phi_near(const double phi, const double reference)
  {
    double out = phi;
    while (out - reference > TMath::Pi()) out -= 2.0 * TMath::Pi();
    while (out - reference < -TMath::Pi()) out += 2.0 * TMath::Pi();
    return out;
  }

  bool is_good_number(const double x)
  {
    return std::isfinite(x) && std::fabs(x) < 1.0e30;
  }

  int track_color(const unsigned int itrk)
  {
    static const int colors[] = {
      kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
      kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7, kTeal + 3
    };
    return colors[itrk % (sizeof(colors) / sizeof(colors[0]))];
  }

  unsigned long long make_unique_hit_id(const TrkrDefs::hitsetkey hsk,
                                        const TrkrDefs::hitkey hk)
  {
    return (static_cast<unsigned long long>(hsk) << 32) |
           static_cast<unsigned long long>(hk);
  }

  void style_fit_line_3d(TPolyLine3D* line, const int color)
  {
    if (!line) return;
    line->SetLineColor(color);
    line->SetLineWidth(4);
    line->SetLineStyle(1);
  }

  double sagitta_model_derivative(const double xrot,
                                  const double x0,
                                  const double invR)
  {
    const double dx = xrot - x0;
    const double dx2 = dx * dx;
    const double invR2 = invR * invR;
    const double invR3 = invR2 * invR;
    const double invR5 = invR3 * invR2;
    return -invR * dx - 0.5 * invR3 * dx2 * dx - 0.375 * invR5 * dx2 * dx2 * dx;
  }

  double sagitta_phi_at_radius(const double radius,
                               const Tpc_FittingTools::SagittaFit& fit)
  {
    const double c = std::cos(fit.theta);
    const double s = std::sin(fit.theta);
    double yy = std::tan(fit.theta) * radius;

    for (unsigned int iter = 0; iter < 25; ++iter)
    {
      const double xrot = c * radius + s * yy;
      const double yrot = -s * radius + c * yy;
      const double f = Tpc_FittingTools::sagittaModel(xrot, fit.S, fit.x0, fit.invR);
      const double g = yrot - f;
      const double df = sagitta_model_derivative(xrot, fit.x0, fit.invR);
      const double dg = c - df * s;
      if (std::fabs(dg) < 1.0e-12) break;
      const double step = g / dg;
      yy -= step;
      if (std::fabs(step) < 1.0e-10) break;
    }

    return fit.b + yy;
  }

  struct RadiusSort
  {
    const std::vector<Tpc_AssembledTrackDisplay::HitPoint>* pts;
    explicit RadiusSort(const std::vector<Tpc_AssembledTrackDisplay::HitPoint>* p) : pts(p) {}
    bool operator()(unsigned int a, unsigned int b) const
    {
      const Tpc_AssembledTrackDisplay::HitPoint& pa = (*pts)[a];
      const Tpc_AssembledTrackDisplay::HitPoint& pb = (*pts)[b];
      if (pa.radius != pb.radius) return pa.radius < pb.radius;
      return pa.tbin < pb.tbin;
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
    double rmin;
    double rmax;
  };

  FitResult fit_assembled_track_points(const std::vector<Tpc_AssembledTrackDisplay::HitPoint>& pts,
                                  const int fit_mode,
                                  const double weight_power,
                                  const double weight_floor_frac)
  {
    FitResult result;
    if (pts.size() < 2) return result;

    std::vector<unsigned int> order;
    order.reserve(pts.size());
    for (unsigned int i = 0; i < pts.size(); ++i) order.push_back(i);
    std::sort(order.begin(), order.end(), RadiusSort(&pts));

    double max_adc = 0.0;
    for (const auto& p : pts) max_adc = std::max(max_adc, static_cast<double>(p.adc));
    if (max_adc <= 0.0) max_adc = 1.0;

    std::vector<Tpc_FittingTools::FitPoint> radius_phi_points;
    std::vector<Tpc_FittingTools::FitPoint> radius_tbin_points;
    radius_phi_points.reserve(pts.size());
    radius_tbin_points.reserve(pts.size());

    bool first = true;
    double phi_reference = 0.0;
    result.rmin = std::numeric_limits<double>::max();
    result.rmax = -std::numeric_limits<double>::max();

    for (unsigned int io = 0; io < order.size(); ++io)
    {
      const Tpc_AssembledTrackDisplay::HitPoint& p = pts[order[io]];
      if (!is_good_number(p.radius) || !is_good_number(p.global_phi)) continue;

      double phi = p.global_phi;
      if (first)
      {
        phi_reference = phi;
        first = false;
      }
      else
      {
        phi = unwrap_phi_near(phi, phi_reference);
        phi_reference = phi;
      }

      const double w = Tpc_FittingTools::adcWeight(static_cast<double>(p.adc), max_adc,
                                         weight_power, weight_floor_frac);
      radius_phi_points.emplace_back(p.radius, phi, w);
      radius_tbin_points.emplace_back(p.radius, static_cast<double>(p.tbin), w);
      result.rmin = std::min(result.rmin, p.radius);
      result.rmax = std::max(result.rmax, p.radius);
    }

    if (radius_phi_points.size() < 2 || radius_tbin_points.size() < 2) return result;

    const Tpc_FittingTools::LineFit phi_fit = Tpc_FittingTools::fitLine(radius_phi_points);
    const Tpc_FittingTools::LineFit tbin_fit = Tpc_FittingTools::fitLine(radius_tbin_points);
    if (!phi_fit.ok || !tbin_fit.ok) return result;

    result.ok = true;
    result.phi_slope = phi_fit.slope;
    result.phi_intercept = phi_fit.intercept;
    result.tbin_slope = tbin_fit.slope;
    result.tbin_intercept = tbin_fit.intercept;

    if (fit_mode == Tpc_FittingTools::FIT_SAGITTA && radius_phi_points.size() >= 3)
    {
      result.phi_sagitta = Tpc_FittingTools::fitSagitta(radius_phi_points);
      result.use_sagitta = result.phi_sagitta.ok;
    }

    return result;
  }

  double fit_direction(const FitResult& fit)
  {
    if (!fit.ok) return 0.0;
    return fit.use_sagitta ? fit.phi_sagitta.theta : std::atan(fit.phi_slope);
  }

  double fit_theta(const FitResult& fit)
  {
    if (!fit.ok) return 0.0;
    return std::atan(fit.tbin_slope);
  }

  double fit_curvature(const FitResult& fit)
  {
    if (!fit.ok) return 0.0;
    return fit.use_sagitta ? fit.phi_sagitta.invR : 0.0;
  }

  void flush_fit_segment(std::vector<TPolyLine3D*>& lines,
                         std::vector<double>& tbin_values,
                         std::vector<double>& phi_values,
                         std::vector<double>& radius_values,
                         const int color,
                         const bool draw_xy)
  {
    if (tbin_values.size() < 2)
    {
      tbin_values.clear();
      phi_values.clear();
      radius_values.clear();
      return;
    }

    TPolyLine3D* line = new TPolyLine3D(static_cast<int>(tbin_values.size()));
    for (unsigned int i = 0; i < tbin_values.size(); ++i)
    {
      if (draw_xy)
      {
        const double x = radius_values[i] * std::cos(phi_values[i]);
        const double y = radius_values[i] * std::sin(phi_values[i]);
        line->SetPoint(static_cast<int>(i), tbin_values[i], x, y);
      }
      else
      {
        line->SetPoint(static_cast<int>(i), tbin_values[i], phi_values[i], radius_values[i]);
      }
    }
    style_fit_line_3d(line, color);
    lines.push_back(line);

    tbin_values.clear();
    phi_values.clear();
    radius_values.clear();
  }

  void add_fit_lines(std::vector<TPolyLine3D*>& lines,
                     const FitResult& fit,
                     const int color,
                     const bool draw_xy)
  {
    if (!fit.ok || fit.rmax <= fit.rmin) return;

    const int npts = 120;
    std::vector<double> seg_tbin;
    std::vector<double> seg_phi;
    std::vector<double> seg_radius;
    double previous_phi = 0.0;
    bool have_previous = false;

    for (int i = 0; i < npts; ++i)
    {
      const double f = static_cast<double>(i) / static_cast<double>(npts - 1);
      const double radius = fit.rmin + f * (fit.rmax - fit.rmin);
      double phi = fit.use_sagitta ? sagitta_phi_at_radius(radius, fit.phi_sagitta)
                                   : fit.phi_slope * radius + fit.phi_intercept;
      phi = wrap_phi(phi);
      const double tbin = fit.tbin_slope * radius + fit.tbin_intercept;

      if (!is_good_number(phi) || !is_good_number(tbin) || !is_good_number(radius)) continue;

      const bool crosses_phi_edge = have_previous && std::fabs(phi - previous_phi) > TMath::Pi();
      if (crosses_phi_edge)
      {
        flush_fit_segment(lines, seg_tbin, seg_phi, seg_radius, color, draw_xy);
      }

      seg_tbin.push_back(tbin);
      seg_phi.push_back(phi);
      seg_radius.push_back(radius);
      previous_phi = phi;
      have_previous = true;
    }

    flush_fit_segment(lines, seg_tbin, seg_phi, seg_radius, color, draw_xy);
  }
}

Tpc_AssembledTrackDisplay::HitPoint::HitPoint()
  : ok(false)
  , hitsetkey(0)
  , hitkey(0)
  , layer(0)
  , region(0)
  , sector(0)
  , side(0)
  , pad(0)
  , tbin(0)
  , adc(0)
  , global_phi(0.0)
  , radius(0.0)
{
}

Tpc_AssembledTrackDisplay::Tpc_AssembledTrackDisplay(const std::string& name,
                                   const std::string& outfilename,
                                   const std::string& trackNodeName,
                                   const unsigned int maxEventDisplays)
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
  , m_plotMinModules(0)
  , m_plotMaxModules(999999)
  , m_plotMinDirection(-1.0e30)
  , m_plotMaxDirection(1.0e30)
  , m_plotMinTheta(-1.0e30)
  , m_plotMaxTheta(1.0e30)
  , m_plotMinCurvature(-1.0e30)
  , m_plotMaxCurvature(1.0e30)
{
}

Tpc_AssembledTrackDisplay::~Tpc_AssembledTrackDisplay()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
}

int Tpc_AssembledTrackDisplay::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << "Tpc_AssembledTrackDisplay::Init - cannot open output file "
              << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_idealPadMap) m_idealPadMap = new IdealPadMap();
  if (!m_idealPadMap->is_loaded() && m_idealPadMap->load_from_cdb(Verbosity()) != 0)
  {
    std::cerr << "Tpc_AssembledTrackDisplay::Init - failed to load IdealPadMap from CDB" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_outfile->mkdir("events");
  std::cout << "Tpc_AssembledTrackDisplay::Init - writing up to "
            << m_maxEventDisplays << " events to "
            << m_outfilename << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_evt;

  if (!get_nodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;
  if (!m_outfile || m_eventsSaved >= m_maxEventDisplays) return Fun4AllReturnCodes::EVENT_OK;

  const unsigned int ntracks = m_tracks ? m_tracks->size() : 0;

  TDirectory* eventsTop = m_outfile->GetDirectory("events");
  if (!eventsTop) eventsTop = m_outfile->mkdir("events");
  eventsTop->cd();

  TDirectory* eventDir = eventsTop->mkdir(Form("event_%06u", m_evt));
  if (!eventDir)
  {
    std::cerr << "Tpc_AssembledTrackDisplay::process_event - failed to create event directory" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  eventDir->cd();

  TH3D* h3[2] = {nullptr, nullptr};
  TH3D* h3xy[2] = {nullptr, nullptr};
  TH3D* h3_single[2] = {nullptr, nullptr};
  TH3D* h3xy_single[2] = {nullptr, nullptr};

  std::vector<TPolyLine3D*> fit_lines_tpr[2];
  std::vector<TPolyLine3D*> fit_lines_txy[2];
  std::vector<TPolyLine3D*> fit_lines_tpr_single[2];
  std::vector<TPolyLine3D*> fit_lines_txy_single[2];

  std::set<unsigned long long> filled_hit_ids[2];
  std::set<unsigned long long> filled_hit_ids_single[2];

  for (unsigned int side = 0; side < 2; ++side)
  {
    h3[side] = new TH3D(Form("h3_evt%06u_tpc_assembledtrack_hits_side%u", m_evt, side),
                        Form("event %u side %u assembled tracks;timebin;global #phi;radius [cm]", m_evt, side),
                        512, -0.5, 511.5,
                        720, -TMath::Pi(), TMath::Pi(),
                        100, 30, 80.0);
    h3[side]->SetStats(0);
    h3[side]->SetDirectory(nullptr);

    h3xy[side] = new TH3D(Form("h3_evt%06u_tpc_assembledtrack_hits_xy_side%u", m_evt, side),
                          Form("event %u side %u assembled tracks;timebin;x [cm];y [cm]", m_evt, side),
                          512, -0.5, 511.5,
                          160, -80.0, 80.0,
                          160, -80.0, 80.0);
    h3xy[side]->SetStats(0);
    h3xy[side]->SetDirectory(nullptr);

    h3_single[side] = new TH3D(Form("h3_evt%06u_single_module_tpc_assembledtrack_hits_side%u", m_evt, side),
                               Form("event %u side %u single-module assembled tracks;timebin;global #phi;radius [cm]", m_evt, side),
                               512, -0.5, 511.5,
                               720, -TMath::Pi(), TMath::Pi(),
                               100, 30, 80.0);
    h3_single[side]->SetStats(0);
    h3_single[side]->SetDirectory(nullptr);

    h3xy_single[side] = new TH3D(Form("h3_evt%06u_single_module_tpc_assembledtrack_hits_xy_side%u", m_evt, side),
                                 Form("event %u side %u single-module assembled tracks;timebin;x [cm];y [cm]", m_evt, side),
                                 512, -0.5, 511.5,
                                 160, -80.0, 80.0,
                                 160, -80.0, 80.0);
    h3xy_single[side]->SetStats(0);
    h3xy_single[side]->SetDirectory(nullptr);
  }

  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const Tpc_AssembledTrack* trk = m_tracks->get_track(itrk);
    if (!trk) continue;

    const int trk_side = trk->get_side();
    if (trk_side < 0 || trk_side > 1) continue;
    const unsigned int side = static_cast<unsigned int>(trk_side);
    const bool single_module_track = (trk->get_nsegments() == 1U);

    const unsigned int nmodules = trk->get_nsegments();
    if (nmodules < m_plotMinModules || nmodules > m_plotMaxModules) continue;

    std::vector<HitPoint> pts;
    pts.reserve(trk->size_hit_indices());

    for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
    {
      const Tpc_AssembledTrack::HitIndex idx = trk->get_hit_index(ih);
      const HitPoint p = make_hit_point(idx.first, idx.second);
      if (!p.ok) continue;
      pts.push_back(p);
    }

    const FitResult fit = fit_assembled_track_points(pts, m_fitMode, m_fitWeightPower, m_fitWeightFloorFrac);
    if (!fit.ok) continue;

    const double direction = fit_direction(fit);
    const double theta = fit_theta(fit);
    const double curvature = fit_curvature(fit);
    if (direction < m_plotMinDirection || direction > m_plotMaxDirection) continue;
    if (theta < m_plotMinTheta || theta > m_plotMaxTheta) continue;
    if (curvature < m_plotMinCurvature || curvature > m_plotMaxCurvature) continue;

    for (unsigned int ip = 0; ip < pts.size(); ++ip)
    {
      const HitPoint& p = pts[ip];
      const unsigned long long uid = make_unique_hit_id(p.hitsetkey, p.hitkey);
      const double phi = wrap_phi(p.global_phi);
      const double x = p.radius * std::cos(p.global_phi);
      const double y = p.radius * std::sin(p.global_phi);

      if (filled_hit_ids[side].insert(uid).second)
      {
        h3[side]->Fill(static_cast<double>(p.tbin), phi, p.radius, static_cast<double>(p.adc));
        h3xy[side]->Fill(static_cast<double>(p.tbin), x, y, static_cast<double>(p.adc));
      }

      if (single_module_track && filled_hit_ids_single[side].insert(uid).second)
      {
        h3_single[side]->Fill(static_cast<double>(p.tbin), phi, p.radius, static_cast<double>(p.adc));
        h3xy_single[side]->Fill(static_cast<double>(p.tbin), x, y, static_cast<double>(p.adc));
      }
    }

    const int color = track_color(itrk);
    add_fit_lines(fit_lines_tpr[side], fit, color, false);
    add_fit_lines(fit_lines_txy[side], fit, color, true);

    if (single_module_track)
    {
      add_fit_lines(fit_lines_tpr_single[side], fit, color, false);
      add_fit_lines(fit_lines_txy_single[side], fit, color, true);
    }
  }

  for (unsigned int side = 0; side < 2; ++side)
  {
    TCanvas* c3 = new TCanvas(Form("c3_evt%06u_timebin_phi_radius_fits_side%u", m_evt, side),
                              Form("event %u side %u assembled-track hits and display fits", m_evt, side),
                              1200, 900);
    h3[side]->Draw("BOX2Z");
    for (unsigned int iline = 0; iline < fit_lines_tpr[side].size(); ++iline)
    {
      if (fit_lines_tpr[side][iline]) fit_lines_tpr[side][iline]->Draw("same");
    }
    c3->Modified();
    c3->Update();
    c3->Write();

    TCanvas* c3xy = new TCanvas(Form("c3_evt%06u_timebin_x_y_fits_side%u", m_evt, side),
                                Form("event %u side %u assembled-track hits and display fits", m_evt, side),
                                1200, 900);
    h3xy[side]->Draw("BOX2Z");
    for (unsigned int iline = 0; iline < fit_lines_txy[side].size(); ++iline)
    {
      if (fit_lines_txy[side][iline]) fit_lines_txy[side][iline]->Draw("same");
    }
    c3xy->Modified();
    c3xy->Update();
    c3xy->Write();

    TCanvas* c3_single = new TCanvas(Form("c3_evt%06u_timebin_phi_radius_fits_single_module_side%u", m_evt, side),
                                     Form("event %u side %u single-module assembled-track hits and display fits", m_evt, side),
                                     1200, 900);
    h3_single[side]->Draw("BOX2Z");
    for (unsigned int iline = 0; iline < fit_lines_tpr_single[side].size(); ++iline)
    {
      if (fit_lines_tpr_single[side][iline]) fit_lines_tpr_single[side][iline]->Draw("same");
    }
    c3_single->Modified();
    c3_single->Update();
    c3_single->Write();
/*
    TCanvas* c3xy_single = new TCanvas(Form("c3_evt%06u_timebin_x_y_fits_single_module_side%u", m_evt, side),
                                       Form("event %u side %u single-module assembled-track hits and display fits", m_evt, side),
                                       1200, 900);
    h3xy_single[side]->Draw("BOX2Z");
    for (unsigned int iline = 0; iline < fit_lines_txy_single[side].size(); ++iline)
    {
      if (fit_lines_txy_single[side][iline]) fit_lines_txy_single[side][iline]->Draw("same");
    }
    c3xy_single->Modified();
    c3xy_single->Update();
    c3xy_single->Write();
    */
  }

  std::cout << "Tpc_AssembledTrackDisplay - saved event " << m_evt
            << " with " << ntracks << " assembled tracks" << std::endl;

  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackDisplay::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
  }

  std::cout << "Tpc_AssembledTrackDisplay::End - events seen: " << m_evt
            << ", events written: " << m_eventsSaved
            << ", output file: " << m_outfilename << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

bool Tpc_AssembledTrackDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_tracks = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_trackNodeName.c_str());

  if (!m_tracks)
  {
    const char* candidate_names[] = {
      "TPC_ASSEMBLEDTRACKS",
      "Tpc_AssembledTracks",
      "Tpc_AssembledTrackContainer",
      "TPC_ASSEMBLEDTRACKCONTAINER",
      "TPC_ASSEMBLEDTRACKS_CONTAINER"
    };

    for (unsigned int i = 0;
         i < sizeof(candidate_names) / sizeof(candidate_names[0]) && !m_tracks;
         ++i)
    {
      m_tracks = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, candidate_names[i]);
      if (m_tracks) m_trackNodeName = candidate_names[i];
    }
  }

  if (!m_idealPadMap) m_idealPadMap = new IdealPadMap();

  if (!m_tracks)
  {
    std::cerr << "Tpc_AssembledTrackDisplay - could not find Tpc_AssembledTrackContainer node" << std::endl;
    return false;
  }

  if (!m_hits)
  {
    std::cerr << "Tpc_AssembledTrackDisplay - missing TRKR_HITSET" << std::endl;
    return false;
  }

  if (!m_idealPadMap || !m_idealPadMap->is_loaded())
  {
    std::cerr << "Tpc_AssembledTrackDisplay - IdealPadMap is not loaded" << std::endl;
    return false;
  }

  return true;
}

Tpc_AssembledTrackDisplay::HitPoint
Tpc_AssembledTrackDisplay::make_hit_point(const TrkrDefs::hitsetkey hsk,
                                 const TrkrDefs::hitkey hk) const
{
  HitPoint p;

  TrkrHitSet* hitset = m_hits ? m_hits->findHitSet(hsk) : nullptr;
  if (!hitset) return p;

  TrkrHit* hit = hitset->getHit(hk);
  if (!hit) return p;

  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.layer = TrkrDefs::getLayer(hsk);
  p.side = static_cast<int>(TpcDefs::getSide(hsk));
  p.pad = TpcDefs::getPad(hk);
  p.tbin = TpcDefs::getTBin(hk);
  p.adc = hit->getAdc();

  if (p.layer < 7 || p.layer > 54) return p;
  p.region = static_cast<unsigned int>((p.layer - 7) / 16);

  if (!m_idealPadMap) return p;
  const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(p.layer);
  if (pads_per_sector == 0U) return p;

  p.sector = (p.pad / pads_per_sector) % 12U;
  p.radius = m_idealPadMap->get_radius(p.layer);
  p.global_phi = wrap_phi(m_idealPadMap->get_phi(static_cast<unsigned int>(p.side), p.layer, p.pad));

  if (!is_good_number(p.radius) || !is_good_number(p.global_phi)) return p;

  p.ok = true;
  return p;
}
