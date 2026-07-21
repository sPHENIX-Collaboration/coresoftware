#include "MicromegasDriftEvaluator.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TVector3.h>

#include <cassert>
#include <climits>
#include <cmath>
#include <format>
#include <iostream>
#include <memory>

namespace
{

  template <class T>
  class range_adaptor
  {
   public:
    explicit range_adaptor(const T& range)
      : m_range(range)
    {
    }
    const typename T::first_type& begin() { return m_range.first; }
    const typename T::second_type& end() { return m_range.second; }

   private:
    T m_range;
  };

  template <class T>
  constexpr T square(T x)
  {
    return x * x;
  }
  template <class T>
  inline T get_r(T x, T y)
  {
    return std::sqrt(square(x) + square(y));
  }

  double normalize_angle(double phi)
  {
    while (phi < 0)
    {
      phi += 2 * M_PI;
    }
    while (phi >= 2 * M_PI)
    {
      phi -= 2 * M_PI;
    }
    return phi;
  }

  bool phi_in_range(double phi, double min, double max)
  {
    phi = normalize_angle(phi);
    min = normalize_angle(min);
    max = normalize_angle(max);
    return (min < max) ? (phi >= min && phi <= max)
                       : (phi >= min || phi <= max);
  }

  //  This function is identical to the version in MicromegasTrackEvaluator_hp.cc

  bool helix_plane_intersection(
      double t_min,
      double t_max,
      double zmin,
      double zmax,
      double R,
      double X0,
      double Y0,
      double intersect_rz,
      double slope_rz,
      const TVector3& ptile,
      const TVector3& ntile,
      TVector3& intersect)
  {
    // Number of iterations and tolerance for Newton Raphson method
    const int max_iter = 10;
    const double tol = 1e-6;

    // Define C
    double C = ntile.X() * (X0 - ptile.X()) + ntile.Y() * (Y0 - ptile.Y()) + ntile.Z() * (intersect_rz - ptile.Z());

    // Defines the function and the corresponding derivative to be used in the Newton Raphson method
    auto f = [&](double t)
    {
      double xt = X0 + R * std::cos(t);
      double yt = Y0 + R * std::sin(t);
      double Rt = std::sqrt(xt * xt + yt * yt);
      return ntile.X() * R * std::cos(t) + ntile.Y() * R * std::sin(t) + ntile.Z() * slope_rz * Rt + C;
    };

    auto df = [&](double t)
    {
      double xt = X0 + R * std::cos(t);
      double yt = Y0 + R * std::sin(t);
      double Rt = std::sqrt(xt * xt + yt * yt);
      return -ntile.X() * R * std::sin(t) + ntile.Y() * R * std::cos(t) + ntile.Z() * R * slope_rz * (Y0 * std::cos(t) - X0 * std::sin(t)) / Rt;
    };

    auto solve_from = [&](double t_seed, TVector3& result) -> bool
    {
      double t = t_seed;
      for (int i = 0; i < max_iter; ++i)
      {
        double ft = f(t);
        double dft = df(t);
        if (std::abs(dft) < 1e-8)
        {
          return false;
        }
        double t_new = t - ft / dft;

        double x = X0 + R * std::cos(t_new);
        double y = Y0 + R * std::sin(t_new);
        double Rt_n = std::sqrt(x * x + y * y);
        double z = slope_rz * Rt_n + intersect_rz;
        double phi = std::atan2(y, x);

        TVector3 cand(x, y, z);
        bool phi_ok = phi_in_range(phi, t_min - 1e-4, t_max + 1e-4);
        bool z_ok = (z >= zmin - 1e-4 && z <= zmax + 1e-4);
        bool proj_ok = (std::abs(ntile.Dot(cand - ptile)) <= 0.05);

        if (std::abs(t_new - t) < tol && proj_ok && phi_ok && z_ok)
        {
          result = cand;
          return true;
        }
        t = t_new;
      }
      return false;
    };

    auto wrap = [&](double t)
    {
      while (t > t_max)
      {
        t -= 2 * M_PI;
      }
      while (t < t_min)
      {
        t += 2 * M_PI;
      }
      return t;
    };

    std::vector<double> t_seeds;
    double t_center = 0.5 * (t_min + t_max);
    double delta = 2.0 * M_PI / 3.0;

    // Wrap the angle
    for (int i = 0; i < 3; ++i)
    {
      double t = wrap(t_center + i * delta);
      t_seeds.push_back(t);
    }

    // Looks for the solution within the tile acceptance in three different phi seeds in the Newton-Raphson (helix_plane could have more than one solution)
    for (double t_seed : t_seeds)
    {
      if (solve_from(t_seed, intersect))
      {
        return true;
      }
    }
    return false;
  }

  // this is a piecewise fit function for the drift velocity plot
  double fit_function_2d(double* x, double* par)
  {
    const int itile = static_cast<int>(std::floor(x[0]));
    const double z = x[1];
    if (itile < 0 || itile >= 8)
    {
      TF2::RejectPoint();
      return 0.;
    }
    return par[itile + 1] + par[0] * z;
  }

  double linear_function(double* x, double* par)
  {
    return par[0] * x[0] + par[1];
  }

  const std::array<const char*, 8> k_tile_names =
      {"SCOZ", "SCIZ", "NCIZ", "NCOZ", "SEZ", "NEZ", "SWZ", "NWZ"};

}  // namespace

MicromegasDriftEvaluator::MicromegasDriftEvaluator(const std::string& name)
  : SubsysReco(name)
{
}

// ---------------------------------------------------------------------------
int MicromegasDriftEvaluator::Init(PHCompositeNode* topNode)
{
  std::cout << Name() << "::Init"
            << " drift_velocity=" << m_drift_velocity << " cm/ns"
            << " min_tpc_layer=" << m_min_tpc_layer
            << " max_tpc_layer=" << m_max_tpc_layer
            << std::endl;

  PHNodeIterator iter(topNode);
  auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::Init - DST node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  iter = PHNodeIterator(dstNode);
  auto* evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if (!evalNode)
  {
    evalNode = new PHCompositeNode("EVAL");
    dstNode->addNode(evalNode);
  }

  auto* newNode = new PHIODataNode<PHObject>(new Container, "MicromegasDriftEvaluator::Container", "PHObject");
  newNode->SplitLevel(99);
  evalNode->addNode(newNode);

  m_hist3D = new TH3F("MicromegasDriftEval_hist3D", ";tile;z_{track} (cm);#Deltaz (track#minuscluster) (cm)", 8, 0, 8, 220, -110, 110, 100, -10, 10);
  m_hist3D->SetDirectory(nullptr);

  return Fun4AllReturnCodes::EVENT_OK;
}

// ---------------------------------------------------------------------------
int MicromegasDriftEvaluator::InitRun(PHCompositeNode* topNode)
{
  return load_nodes(topNode);
}

// ---------------------------------------------------------------------------
int MicromegasDriftEvaluator::process_event(PHCompositeNode* topNode)
{
  const auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK)
  {
    return res;
  }

  if (m_container)
  {
    m_container->Reset();
  }
  evaluate_tracks();

  return Fun4AllReturnCodes::EVENT_OK;
}

// ---------------------------------------------------------------------------
int MicromegasDriftEvaluator::End(PHCompositeNode* /*topNode*/)
{
  if (!m_hist3D)
  {
    std::cerr << Name() << "::End - histogram not found, skipping fit." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const int nEntries = static_cast<int>(m_hist3D->GetEntries());
  std::cout << Name() << "::End - fitting " << nEntries << " entries" << std::endl;

  auto* h_fit = new TH2F("h_fit", "", 8, 0, 8, 220, -110, 110);
  h_fit->SetDirectory(nullptr);

  for (int j = 0; j < 8; ++j)
  {
    m_hist3D->GetXaxis()->SetRange(j + 1, j + 1);
    auto* h2d = static_cast<TH2F*>(m_hist3D->Project3D("zy"));
    h2d->SetName(std::format("h_{}", k_tile_names[j]).c_str());
    h2d->SetDirectory(nullptr);

    // Fit vertical slices; require a minimum of 10 entries per slice
    h2d->FitSlicesY(nullptr, 0, -1, 10);
    auto* h_mean = static_cast<TH1F*>(gDirectory->Get(std::format("h_{}_1", k_tile_names[j]).c_str()));

    if (!h_mean)
    {
      delete h2d;
      continue;
    }

    for (int i = 0; i < h_mean->GetNbinsX(); ++i)
    {
      const double entries = h2d->Integral(i + 1, i + 1, 1, m_hist3D->GetNbinsZ());
      if (entries > 0)
      {
        h_fit->SetBinContent(j + 1, i + 1, h_mean->GetBinContent(i + 1));
      }
    }
    delete h2d;
  }

  m_hist3D->GetXaxis()->SetRange(0, 0);

  // This part fits the 8 micromegas modules one at a time but constrains the slopes to be identical. This eliminates the need for perfect translational TPOT alignment
  auto* fit2d = new TF2("fit2d", fit_function_2d, 0, 8, -110, 110, 9);
  for (int i = 0; i < 9; ++i)
  {
    fit2d->SetParameter(i, 0.0);
  }

  h_fit->Fit(fit2d, "0R");

  const double slope = fit2d->GetParameter(0);
  const double slope_err = fit2d->GetParError(0);
  const double new_drift = m_drift_velocity / (1.0 + slope);
  const double drift_err = m_drift_velocity / std::pow(1.0 + slope, 2) * slope_err;

  std::cout << Name() << "::End" << " slope=" << slope << " input_drift=" << m_drift_velocity << " cm/ns" << " new_drift=" << new_drift << " cm/ns +/- " << drift_err << " cm/ns" << std::endl;

  // Plot the whole thing
  auto* canvas = new TCanvas("drift_calib", "Drift velocity calibration", 2000, 1000);
  canvas->Divide(4, 2);

  for (int j = 0; j < 8; ++j)
  {
    canvas->cd(j + 1);

    m_hist3D->GetXaxis()->SetRange(j + 1, j + 1);
    auto* h2d = static_cast<TH2F*>(m_hist3D->Project3D("zy"));
    h2d->SetName(std::format("hplot_{}", k_tile_names[j]).c_str());
    h2d->SetTitle(std::format("{};z_{{track}} (cm);#Deltaz (track#minuscluster) (cm)", k_tile_names[j]).c_str());
    h2d->SetStats(false);
    h2d->Draw("COLZ");

    auto* h_fit_proj = h_fit->ProjectionY(std::format("h_fit_proj_{}", j).c_str(), j + 1, j + 1);  // These give you the Gaussian means for each slice
    h_fit_proj->SetMarkerStyle(20);
    h_fit_proj->SetMarkerColor(kRed);
    h_fit_proj->SetLineColor(kBlack);

    auto* f1d = new TF1(std::format("f1d_{}", j).c_str(), linear_function, -110, 110, 2);
    f1d->SetParameter(0, slope);
    f1d->SetParameter(1, fit2d->GetParameter(j + 1));
    f1d->SetLineColor(kGreen + 2);
    f1d->SetLineWidth(2);
    f1d->Draw("same");

    auto* leg = new TLegend(0.35, 0.75, 0.92, 0.92);
    leg->SetHeader(std::format("{} entries, v_{{in}}={:.2f} m/ms", nEntries, m_drift_velocity * 1e4).c_str(), "C");
    leg->AddEntry(h_fit_proj, "Gaussian slice mean", "p");
    leg->AddEntry(f1d, std::format("slope={:.4f}  v_{{new}}={:.3f}#pm{:.3f} m/ms", slope, new_drift * 1e4, drift_err * 1e4).c_str(), "l");
    leg->Draw();
  }

  m_hist3D->GetXaxis()->SetRange(0, 0);

  canvas->SaveAs(m_plot_filename.c_str());
  std::cout << Name() << "::End - QA canvas saved to " << m_plot_filename << std::endl;

  // write histograms, fit and results to a ROOT file
  if (!m_root_filename.empty())
  {
    std::unique_ptr<TFile> outfile(TFile::Open(m_root_filename.c_str(), "RECREATE"));
    if (outfile && !outfile->IsZombie())
    {
      outfile->cd();
      m_hist3D->Write();
      h_fit->Write("h_fit_micromegas");
      fit2d->Write();
      canvas->Write();
      TParameter<double>("slope", slope).Write();
      TParameter<double>("drift_velocity_in", m_drift_velocity).Write();
      TParameter<double>("drift_velocity_new", new_drift).Write();
      TParameter<double>("drift_velocity_err", drift_err).Write();
      outfile->Close();
      std::cout << Name() << "::End - histograms and fit results saved to " << m_root_filename << std::endl;
    }
    else
    {
      std::cerr << Name() << "::End - could not open " << m_root_filename << " for writing." << std::endl;
    }
  }

  delete canvas;
  delete fit2d;
  delete h_fit;

  return Fun4AllReturnCodes::EVENT_OK;
}

// ---------------------------------------------------------------------------
int MicromegasDriftEvaluator::load_nodes(PHCompositeNode* topNode)
{
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert(m_tGeometry);

  m_micromegas_geomcontainer = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  assert(m_micromegas_geomcontainer);

  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert(m_cluster_map);

  m_container = findNode::getClass<Container>(topNode, "MicromegasDriftEvaluator::Container");
  assert(m_container);

  m_globalPositionWrapper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

// ---------------------------------------------------------------------------
void MicromegasDriftEvaluator::evaluate_tracks()
{
  if (!(m_tGeometry && m_micromegas_geomcontainer && m_track_map && m_cluster_map && m_container && m_hist3D))
  {
    return;
  }

  m_container->clear_tracks();

  for (const auto& [track_id, track] : *m_track_map)
  {
    // valid crossing
    const auto crossing = track->get_crossing();
    if (crossing == SHRT_MAX)
    {
      continue;
    }

    std::vector<Acts::Vector3> tpc_positions;

    // Also count clusters per subsystem for the cuts
    unsigned int n_tpc = 0;
    unsigned int n_mvtx = 0;
    unsigned int n_intt = 0;
    unsigned int n_mm = 0;

    for (const auto* seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (!seed)
      {
        continue;
      }
      for (auto it = seed->begin_cluster_keys(); it != seed->end_cluster_keys(); ++it)
      {
        const auto ckey = *it;
        const auto detid = TrkrDefs::getTrkrId(ckey);
        const auto layer = TrkrDefs::getLayer(ckey);

        switch (detid)
        {
        case TrkrDefs::tpcId:
          ++n_tpc;
          if (layer >= m_min_tpc_layer && layer < m_max_tpc_layer)
          {
            auto* const cl = m_cluster_map->findCluster(ckey);
            if (cl)
            {
              tpc_positions.push_back(
                  m_globalPositionWrapper.getGlobalPositionDistortionCorrected(
                      ckey, cl, crossing));
            }
          }
          break;
        case TrkrDefs::mvtxId:
          ++n_mvtx;
          break;
        case TrkrDefs::inttId:
          ++n_intt;
          break;
        case TrkrDefs::micromegasId:
          ++n_mm;
          break;
        default:
          break;
        }
      }
    }

    // need at least 3 TPC clusters in range
    if (tpc_positions.size() < 3)
    {
      continue;
    }

    const auto [slope_rz, intersect_rz] = TrackFitUtils::line_fit(tpc_positions);
    const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin(tpc_positions);

    // reject badly reconstructed / low-pT tracks
    if (R < 40.0)
    {
      continue;
    }

    const auto mm_range = m_micromegas_geomcontainer->get_begin_end();
    for (const auto& [mm_layer, base_layergeom] : range_adaptor(mm_range))
    {
      const auto* layergeom = static_cast<const CylinderGeomMicromegas*>(base_layergeom);
      assert(layergeom);

      // skip the phi layer. Only the z-view layer matters here
      if (layergeom->get_segmentation_type() !=
          MicromegasDefs::SegmentationType::SEGMENTATION_Z)
      {
        continue;
      }

      const double layer_radius = layergeom->get_radius();
      auto [xplus, yplus, xminus, yminus] =
          TrackFitUtils::circle_circle_intersection(layer_radius, R, X0, Y0);

      if (!std::isfinite(xplus))
      {
        continue;
      }

      // pick the solution closest in phi to the last TPC cluster
      const double last_phi = std::atan2(tpc_positions.back().y(), tpc_positions.back().x());
      const double phi_plus = std::atan2(yplus, xplus);
      const double phi_minus = std::atan2(yminus, xminus);
      const double phi = (std::abs(last_phi - phi_plus) < std::abs(last_phi - phi_minus)) ? phi_plus : phi_minus;

      const double r_cyl = layer_radius;
      const double z_cyl = intersect_rz + slope_rz * r_cyl;
      const TVector3 world_cyl(r_cyl * std::cos(phi), r_cyl * std::sin(phi), z_cyl);

      const int tileid = layergeom->find_tile_cylindrical(world_cyl);
      if (tileid < 0)
      {
        continue;
      }

      const auto tile_center = layergeom->get_world_from_local_coords(tileid, m_tGeometry, {0, 0});
      const TVector3 ptile(tile_center.x(), tile_center.y(), tile_center.z());

      const auto tile_norm = layergeom->get_world_from_local_vect(tileid, m_tGeometry, {0, 0, 1});
      const TVector3 ntile(tile_norm.x(), tile_norm.y(), tile_norm.z());

      const auto phi_range = layergeom->get_phi_range(tileid, m_tGeometry);
      const double zmin = layergeom->get_zmin();
      const double zmax = layergeom->get_zmax();

      TVector3 intersection;
      if (!helix_plane_intersection(phi_range.first, phi_range.second, zmin, zmax, R, X0, Y0, intersect_rz, slope_rz, ptile, ntile, intersection))
      {
        continue;
      }

      const auto local_intersection = layergeom->get_local_from_world_coords(tileid, m_tGeometry, {intersection.x(), intersection.y(), intersection.z()});
      const double y_local = local_intersection.y();

      if (std::abs(y_local) > m_y_local_cut)
      {
        continue;
      }

      // find the nearest TPOT cluster
      const auto hitsetkey = MicromegasDefs::genHitSetKey(mm_layer, MicromegasDefs::SegmentationType::SEGMENTATION_Z, tileid);
      const auto clusrange = m_cluster_map->getClusters(hitsetkey);

      double dmin = -1;
      ClusterStruct best_cluster;

      for (const auto& [ckey, cl] : range_adaptor(clusrange))
      {
        const double cl_y_local = cl->getLocalY();
        const double d = std::abs(y_local - cl_y_local);
        if (dmin < 0 || d < dmin)
        {
          dmin = d;
          const auto gpos = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cl, crossing);
          best_cluster._layer = mm_layer;
          best_cluster._tile = tileid;
          best_cluster._z = gpos.z();
        }
      }

      // require cluster within the z search window
      if (dmin < 0 || dmin > m_z_search_win)
      {
        continue;
      }

      // fill track struct and histogram
      TrackStruct track_struct;
      track_struct._chisquare = track->get_chisq();
      track_struct._ndf = track->get_ndf();
      track_struct._nclusters_tpc = n_tpc;
      track_struct._nclusters_mvtx = n_mvtx;
      track_struct._nclusters_intt = n_intt;
      track_struct._nclusters_micromegas = n_mm;

      track_struct._trk_state_z._layer = mm_layer;
      track_struct._trk_state_z._tile = tileid;
      track_struct._trk_state_z._z = intersection.z();
      track_struct._trk_state_z._y_local = y_local;

      track_struct._found_cluster_z = best_cluster;

      const double z_track = track_struct._trk_state_z._z;
      const double z_cluster = track_struct._found_cluster_z._z;
      m_hist3D->Fill(tileid + 0.5, z_track, z_track - z_cluster);

      m_container->add_track(track_struct);
      break;
    }
  }
}