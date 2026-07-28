#include "MicromegasDriftQA.h"

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TDirectory.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>

#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <format>
#include <iostream>
#include <string>
#include <vector>

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

  //! helix-plane intersection via Newton-Raphson
  //  identical to the version in MicromegasTrackEvaluator_hp.cc
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
    // number of iterations and tolerance for Newton-Raphson method
    const int max_iter = 10;
    const double tol = 1e-6;

    // define C
    const double C = ntile.X() * (X0 - ptile.X()) + ntile.Y() * (Y0 - ptile.Y()) + ntile.Z() * (intersect_rz - ptile.Z());

    // define the function and the corresponding derivative used in the Newton-Raphson method
    auto f = [&](double t)
    {
      const double xt = X0 + R * std::cos(t);
      const double yt = Y0 + R * std::sin(t);
      const double Rt = std::sqrt(xt * xt + yt * yt);
      return ntile.X() * R * std::cos(t) + ntile.Y() * R * std::sin(t) + ntile.Z() * slope_rz * Rt + C;
    };

    auto df = [&](double t)
    {
      const double xt = X0 + R * std::cos(t);
      const double yt = Y0 + R * std::sin(t);
      const double Rt = std::sqrt(xt * xt + yt * yt);
      return -ntile.X() * R * std::sin(t) + ntile.Y() * R * std::cos(t) + ntile.Z() * R * slope_rz * (Y0 * std::cos(t) - X0 * std::sin(t)) / Rt;
    };

    auto solve_from = [&](double t_seed, TVector3& result) -> bool
    {
      double t = t_seed;
      for (int i = 0; i < max_iter; ++i)
      {
        const double ft = f(t);
        const double dft = df(t);
        if (std::abs(dft) < 1e-8)
        {
          return false;
        }
        const double t_new = t - ft / dft;

        const double x = X0 + R * std::cos(t_new);
        const double y = Y0 + R * std::sin(t_new);
        const double Rt_n = std::sqrt(x * x + y * y);
        const double z = slope_rz * Rt_n + intersect_rz;
        const double phi = std::atan2(y, x);

        const TVector3 cand(x, y, z);
        const bool phi_ok = phi_in_range(phi, t_min - 1e-4, t_max + 1e-4);
        const bool z_ok = (z >= zmin - 1e-4 && z <= zmax + 1e-4);
        const bool proj_ok = (std::abs(ntile.Dot(cand - ptile)) <= 0.05);

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

    // the helix-plane equation can have more than one solution:
    // look for a solution within the tile acceptance from three different phi seeds
    std::vector<double> t_seeds(3);
    const double t_center = 0.5 * (t_min + t_max);
    const double delta = 2.0 * M_PI / 3.0;
    for (int i = 0; i < 3; ++i)
    {
      t_seeds[i]=wrap(t_center + i * delta);
    }

    for (const double t_seed : t_seeds)
    {
      if (solve_from(t_seed, intersect))
      {
        return true;
      }
    }
    return false;
  }

  //! piecewise fit function used for the drift velocity extraction
  //  par[0] = constrained slope, par[1..8] = per-tile offsets
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

  //! z-view tile names
  const std::array<const char*, 8> k_tile_names =
      {"SCOZ", "SCIZ", "NCIZ", "NCOZ", "SEZ", "NEZ", "SWZ", "NWZ"};

  //! number of z bins of the dz vs z histograms
  constexpr int k_nzbins = 220;

  //! z_track range (cm)
  constexpr double k_max_z = 110;

  //! dz range (cm)
  constexpr double k_max_dz = 10;

}  // namespace

//____________________________________________________________________________..
MicromegasDriftQA::MicromegasDriftQA(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MicromegasDriftQA::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity())
  {
    std::cout << Name() << "::InitRun"
              << " drift_velocity=" << m_drift_velocity << " cm/ns"
              << " min_tpc_layer=" << m_min_tpc_layer
              << " max_tpc_layer=" << m_max_tpc_layer
              << std::endl;
  }

  const auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK)
  {
    return res;
  }

  createHistos();

  // reference histograms initialized in header file to histos in HistoManager
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (int itile = 0; itile < 8; itile++)
  {
    h_ztrk_dz[itile] = dynamic_cast<TH2*>(hm->getHisto(std::format("{}ztrk_dz_{}", getHistoPrefix(), k_tile_names[itile])));
  }
  h_dz = dynamic_cast<TH1*>(hm->getHisto(std::format("{}dz", getHistoPrefix())));
  h_tile = dynamic_cast<TH1*>(hm->getHisto(std::format("{}tile", getHistoPrefix())));
  h_ylocal = dynamic_cast<TH1*>(hm->getHisto(std::format("{}ylocal", getHistoPrefix())));
  h_ntracks = dynamic_cast<TH1*>(hm->getHisto(std::format("{}ntracks", getHistoPrefix())));
  h_driftSummary = dynamic_cast<TH1*>(hm->getHisto(std::format("{}driftSummary", getHistoPrefix())));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasDriftQA::process_event(PHCompositeNode* topNode)
{
  const auto res = load_nodes(topNode);
  if (res != Fun4AllReturnCodes::EVENT_OK)
  {
    return res;
  }

  int nmatched = 0;

  for (const auto& [track_id, track] : *m_track_map)
  {
    // require valid beam-crossing
    const auto crossing = track->get_crossing();
    if (crossing == SHRT_MAX)
    {
      continue;
    }

    // collect distortion-corrected TPC cluster positions in the selected layer range
    std::vector<Acts::Vector3> tpc_positions;
    for (const auto* seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (!seed)
      {
        continue;
      }
      for (auto it = seed->begin_cluster_keys(); it != seed->end_cluster_keys(); ++it)
      {
        const auto ckey = *it;
        if (TrkrDefs::getTrkrId(ckey) != TrkrDefs::tpcId)
        {
          continue;
        }
        const auto layer = TrkrDefs::getLayer(ckey);
        if (layer < m_min_tpc_layer || layer >= m_max_tpc_layer)
        {
          continue;
        }
        auto* cl = m_cluster_map->findCluster(ckey);
        if (cl)
        {
          tpc_positions.push_back(
              m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cl, crossing));
        }
      }
    }

    // need at least 3 TPC clusters in range
    if (tpc_positions.size() < 3)
    {
      continue;
    }

    // helix fit: straight line in r-z, circle in x-y
    const auto [slope_rz, intersect_rz] = TrackFitUtils::line_fit(tpc_positions);
    const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin(tpc_positions);

    // reject badly reconstructed / low-pT tracks
    if (R < 40.0)
    {
      continue;
    }

    // extrapolate to the TPOT z-view modules
    const auto mm_range = m_micromegas_geomcontainer->get_begin_end();
    for (const auto& [mm_layer, base_layergeom] : range_adaptor(mm_range))
    {
      const auto* layergeom = static_cast<const CylinderGeomMicromegas*>(base_layergeom);
      assert(layergeom);

      // skip the phi layer; only the z-view layer matters here
      if (layergeom->get_segmentation_type() != MicromegasDefs::SegmentationType::SEGMENTATION_Z)
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
      if (!helix_plane_intersection(phi_range.first, phi_range.second, zmin, zmax,
                                    R, X0, Y0, intersect_rz, slope_rz, ptile, ntile, intersection))
      {
        continue;
      }

      const auto local_intersection = layergeom->get_local_from_world_coords(
          tileid, m_tGeometry, {intersection.x(), intersection.y(), intersection.z()});
      const double y_local = local_intersection.y();

      // reject track states near the tile edge
      if (std::abs(y_local) > m_y_local_cut)
      {
        continue;
      }

      // find the nearest TPOT cluster on this tile
      const auto hitsetkey = MicromegasDefs::genHitSetKey(mm_layer, MicromegasDefs::SegmentationType::SEGMENTATION_Z, tileid);
      const auto clusrange = m_cluster_map->getClusters(hitsetkey);

      double dmin = -1;
      double z_cluster = 0;
      for (const auto& [ckey, cl] : range_adaptor(clusrange))
      {
        const double cl_y_local = cl->getLocalY();
        const double d = std::abs(y_local - cl_y_local);
        if (dmin < 0 || d < dmin)
        {
          dmin = d;
          const auto gpos = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cl, crossing);
          z_cluster = gpos.z();
        }
      }

      // require cluster within the z search window
      if (dmin < 0 || dmin > m_z_search_win)
      {
        continue;
      }

      // fill histograms
      const double z_track = intersection.z();
      const double dz = z_track - z_cluster;

      h_ztrk_dz[tileid]->Fill(z_track, dz);
      h_dz->Fill(dz);
      h_tile->Fill(tileid);
      h_ylocal->Fill(y_local);

      ++nmatched;
      break;
    }
  }

  h_ntracks->Fill(nmatched);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasDriftQA::End(PHCompositeNode* /*topNode*/)
{
  if (!(h_ztrk_dz[0] && h_driftSummary))
  {
    std::cout << PHWHERE << " histograms not found, skipping drift velocity fit." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int nEntries = 0;
  for (const auto* h : h_ztrk_dz)
  {
    nEntries += static_cast<int>(h->GetEntries());
  }
  if (Verbosity())
  {
    std::cout << Name() << "::End - fitting " << nEntries << " entries" << std::endl;
  }

  // record input drift velocity even if the fit is skipped
  h_driftSummary->SetBinContent(3, m_drift_velocity);

  if (nEntries < 8 * m_min_slice_entries)
  {
    std::cout << Name() << "::End - not enough entries (" << nEntries << "), skipping drift velocity fit." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // build mean-dz TH2 via FitSlicesY, one tile at a time
  // x = tile [0,8), y = z_track (cm), content = mean dz (cm)
  auto* h_fit = new TH2F("h_fit_micromegas", "", 8, 0, 8, k_nzbins, -k_max_z, k_max_z);
  h_fit->SetDirectory(nullptr);

  for (int itile = 0; itile < 8; ++itile)
  {
    auto* h2d = h_ztrk_dz[itile];

    // fit vertical slices; require a minimum of m_min_slice_entries per slice
    TObjArray slices;
    slices.SetOwner(kTRUE);
    h2d->FitSlicesY(nullptr, 0, -1, m_min_slice_entries, "QNR", &slices);
    auto* h_mean = dynamic_cast<TH1*>(slices.At(1));
    if (!h_mean)
    {
      continue;
    }

    for (int iz = 1; iz <= h_mean->GetNbinsX(); ++iz)
    {
      const double entries = h2d->Integral(iz, iz, 1, h2d->GetNbinsY());
      if (entries > 0)
      {
        h_fit->SetBinContent(itile + 1, iz, h_mean->GetBinContent(iz));
      }
    }
  }

  // 2D piecewise fit: the eight tiles are fitted simultaneously with a shared
  // slope and per-tile offsets. This eliminates the need for perfect
  // translational TPOT alignment.
  auto* fit2d = new TF2("fit2d_micromegas", fit_function_2d, 0, 8, -k_max_z, k_max_z, 9);
  for (int i = 0; i < 9; ++i)
  {
    fit2d->SetParameter(i, 0.0);
  }
  h_fit->Fit(fit2d, "0RQ");

  const double slope = fit2d->GetParameter(0);
  const double slope_err = fit2d->GetParError(0);
  const double new_drift = m_drift_velocity / (1.0 + slope);
  const double drift_err = m_drift_velocity / std::pow(1.0 + slope, 2) * slope_err;

  std::cout << Name() << "::End"
            << " slope=" << slope
            << " input_drift=" << m_drift_velocity << " cm/ns"
            << " new_drift=" << new_drift << " cm/ns +/- " << drift_err << " cm/ns"
            << std::endl;

  // store fit results in the summary histogram
  h_driftSummary->SetBinContent(1, slope);
  h_driftSummary->SetBinContent(2, slope_err);
  h_driftSummary->SetBinContent(4, new_drift);
  h_driftSummary->SetBinContent(5, drift_err);

  delete fit2d;
  delete h_fit;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MicromegasDriftQA::load_nodes(PHCompositeNode* topNode)
{
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << " ActsGeometry node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_micromegas_geomcontainer = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!m_micromegas_geomcontainer)
  {
    std::cout << PHWHERE << " CYLINDERGEOM_MICROMEGAS_FULL node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);
  if (!m_track_map)
  {
    std::cout << PHWHERE << " " << m_trackmapname << " node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_cluster_map)
  {
    std::cout << PHWHERE << " TRKR_CLUSTER node missing, abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_globalPositionWrapper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
std::string MicromegasDriftQA::getHistoPrefix() const
{
  // define prefix to all histos in HistoManager
  return std::string("h_") + Name() + std::string("_");
}

//____________________________________________________________________________..
void MicromegasDriftQA::createHistos()
{
  // initialize HistoManager
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create and register histos in HistoManager
  for (int itile = 0; itile < 8; itile++)
  {
    auto* h = new TH2F(std::format("{}ztrk_dz_{}", getHistoPrefix(), k_tile_names[itile]).c_str(),
                       std::format("{};z_{{track}} (cm);#Deltaz (track#minuscluster) (cm)", k_tile_names[itile]).c_str(),
                       k_nzbins, -k_max_z, k_max_z, 100, -k_max_dz, k_max_dz);
    hm->registerHisto(h);
  }

  {
    auto* h = new TH1F(std::format("{}dz", getHistoPrefix()).c_str(),
                       ";#Deltaz (track#minuscluster) (cm);track states", 100, -k_max_dz, k_max_dz);
    hm->registerHisto(h);
  }

  {
    auto* h = new TH1F(std::format("{}tile", getHistoPrefix()).c_str(),
                       ";tile;track states", 8, -0.5, 7.5);
    for (int itile = 0; itile < 8; itile++)
    {
      h->GetXaxis()->SetBinLabel(itile + 1, k_tile_names[itile]);
    }
    hm->registerHisto(h);
  }

  {
    auto* h = new TH1F(std::format("{}ylocal", getHistoPrefix()).c_str(),
                       ";y_{local} (cm);track states", 100, -30, 30);
    hm->registerHisto(h);
  }

  {
    auto* h = new TH1F(std::format("{}ntracks", getHistoPrefix()).c_str(),
                       ";matched track states per event;events", 20, -0.5, 19.5);
    hm->registerHisto(h);
  }

  {
    // summary of the drift velocity fit performed in End()
    auto* h = new TH1F(std::format("{}driftSummary", getHistoPrefix()).c_str(),
                       "drift velocity fit summary", 5, 0.5, 5.5);
    h->GetXaxis()->SetBinLabel(1, "slope");
    h->GetXaxis()->SetBinLabel(2, "slope_err");
    h->GetXaxis()->SetBinLabel(3, "v_{in} (cm/ns)");
    h->GetXaxis()->SetBinLabel(4, "v_{new} (cm/ns)");
    h->GetXaxis()->SetBinLabel(5, "v_{new} err (cm/ns)");
    hm->registerHisto(h);
  }
}