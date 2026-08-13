#include "TpcV0CandidateTree.h"

#include <tpctrackreco/TpcTrackHelixFitter.h>
#include <tpctrackreco/TpcTrackKalmanFitter.h>
#include <tpctrackreco/Tpc_PolyCluster.h>
#include <tpctrackreco/Tpc_PolyClusterContainer.h>
#include <tpctrackreco/Tpc_PolyTrack.h>
#include <tpctrackreco/Tpc_PolyTrackContainer.h>
#include <tpctrackreco/Tpc_PolyTrackVertexContainer.h>

#include <trackbase/TrkrDefs.h>

#include <g4main/PHG4EventHeader.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <ffaobjects/EventHeader.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <utility>

namespace
{
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kPionMass = 0.13957039;
  constexpr double kProtonMass = 0.938272088;

  template <class T>
  constexpr T square(const T &value)
  {
    return value * value;
  }

  [[maybe_unused]] int sign_to_charge(const double value)
  {
    if (!std::isfinite(value) || value == 0.0)
    {
      return 0;
    }
    return value > 0.0 ? 1 : -1;
  }

  double normalize_phi(const double phi)
  {
    return std::atan2(std::sin(phi), std::cos(phi));
  }

  double unwrap_to_near(double theta, const double reference)
  {
    while (theta - reference > kPi)
    {
      theta -= 2.0 * kPi;
    }
    while (theta - reference < -kPi)
    {
      theta += 2.0 * kPi;
    }
    return theta;
  }

  bool helix_point_at_beam_radius(const TpcTrackHelix &helix,
                                  const double beam_radius,
                                  const double theta_hint,
                                  double &theta,
                                  TpcTrackVec3 &position)
  {
    if (beam_radius <= 0.0 || helix.radius <= 0.0 ||
        !std::isfinite(beam_radius) || !std::isfinite(helix.radius))
    {
      return false;
    }

    const double center_radius = std::hypot(helix.cx, helix.cy);
    if (center_radius <= 1.0e-12)
    {
      theta = theta_hint;
      position = TpcTrackHelixFitter::point(helix, theta);
      return TpcTrackHelixFitter::finite(position);
    }

    const double rhs = (square(beam_radius) - square(center_radius) - square(helix.radius)) /
                       (2.0 * helix.radius * center_radius);
    if (rhs < -1.0 - 1.0e-9 || rhs > 1.0 + 1.0e-9)
    {
      return false;
    }

    const double clamped_rhs = std::clamp(rhs, -1.0, 1.0);
    const double center_phi = std::atan2(helix.cy, helix.cx);
    const double delta = std::acos(clamped_rhs);
    const double theta_a = unwrap_to_near(center_phi + delta, theta_hint);
    const double theta_b = unwrap_to_near(center_phi - delta, theta_hint);
    theta = (std::abs(theta_a - theta_hint) <= std::abs(theta_b - theta_hint)) ? theta_a : theta_b;
    position = TpcTrackHelixFitter::point(helix, theta);
    return TpcTrackHelixFitter::finite(position);
  }

  bool helix_cluster_residual(const TpcTrackHelix &helix,
                              const TpcTrackPoint &point,
                              double &previous_theta,
                              bool &have_previous_theta,
                              TpcTrackVec3 &fit_position,
                              double &residual_r,
                              double &residual_rphi,
                              double &residual_z)
  {
    const TpcTrackVec3 &cluster = point.position;
    const double cluster_r = std::hypot(cluster.x, cluster.y);
    double theta_hint = std::atan2(cluster.y - helix.cy,
                                   cluster.x - helix.cx);
    theta_hint = unwrap_to_near(theta_hint,
                                have_previous_theta ? previous_theta
                                                    : helix.theta_first);

    double theta = theta_hint;
    if (!helix_point_at_beam_radius(helix, cluster_r, theta_hint, theta, fit_position))
    {
      fit_position = TpcTrackHelixFitter::point(helix, theta_hint);
      theta = theta_hint;
      if (!TpcTrackHelixFitter::finite(fit_position))
      {
        return false;
      }
    }

    previous_theta = theta;
    have_previous_theta = true;

    const double fit_r = std::hypot(fit_position.x, fit_position.y);
    const double cluster_phi = std::atan2(cluster.y, cluster.x);
    const double fit_phi = std::atan2(fit_position.y, fit_position.x);
    residual_r = cluster_r - fit_r;
    residual_rphi = cluster_r * normalize_phi(cluster_phi - fit_phi);
    residual_z = cluster.z - fit_position.z;
    return std::isfinite(residual_r) &&
           std::isfinite(residual_rphi) &&
           std::isfinite(residual_z);
  }
}  // namespace

TpcV0CandidateTree::TpcV0CandidateTree(const std::string &name,
                                       const std::string &filename)
  : SubsysReco(name)
  , m_filename(filename)
{
  m_kalman_config.bfield_t = m_bfield_t;
  m_kalman_config.point_order = m_point_order;
}

void TpcV0CandidateTree::set_primary_vertex(const double x, const double y, const double z)
{
  m_fixed_primary_vertex = {x, y, z};
}

bool TpcV0CandidateTree::set_point_order(const std::string &mode)
{
  PointOrder order = PointOrder::Path;
  if (!parse_point_order(mode, order))
  {
    std::cerr << PHWHERE << Name() << ": unknown point order mode '" << mode
              << "'. Valid modes are path, input, radius, theta-z, auto." << std::endl;
    return false;
  }

  m_point_order = order;
  m_kalman_config.point_order = order;
  return true;
}

bool TpcV0CandidateTree::set_track_fit_method(const std::string &mode)
{
  std::string lowered = mode;
  std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                 [](const unsigned char ch)
                 { return static_cast<char>(std::tolower(ch)); });

  if (lowered == "helix" || lowered == "circle")
  {
    set_fit_helix(true);
    return true;
  }
  if (lowered == "kalman" || lowered == "kf")
  {
    set_fit_kalman(true);
    return true;
  }
  if (lowered == "none" || lowered == "line")
  {
    m_fit_helix_tracks = false;
    m_fit_kalman_tracks = false;
    return true;
  }

  std::cerr << PHWHERE << Name() << ": unknown track fit method '" << mode
            << "'. Valid methods are helix, kalman, and line." << std::endl;
  return false;
}

int TpcV0CandidateTree::Init(PHCompositeNode *topNode)
{
  if (m_use_kalman_field_map && m_kalman_config.magnetic_field == nullptr && topNode != nullptr)
  {
    m_kalman_config.magnetic_field =
        findNode::getClass<PHField>(topNode, PHFieldUtility::GetDSTFieldMapNodeName());
  }

  m_file = new TFile(m_filename.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cout << Name() << ": failed to create output file " << m_filename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_pair_tree = new TTree("pairTree", "TPC truth-point V0 candidates");
  m_track_tree = new TTree("trackTree", "TPC track QA before V0 pairing preselection");
  if (m_write_cluster_residual_tree)
  {
    m_cluster_residual_tree = new TTree("clusterResidualTree", "TPC per-cluster residual QA");
  }
  create_branches();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcV0CandidateTree::process_event(PHCompositeNode *topNode)
{
  const auto event_start = std::chrono::steady_clock::now();
  const double kalman_fit_before = m_timing_kalman_fit_seconds;
  const double kalman_pca_before = m_timing_kalman_pca_seconds;
  const std::uint64_t rkn_propagations_before = m_timing_rkn_propagations;
  const std::uint64_t rkn_steps_before = m_timing_rkn_accepted_steps;
  const std::uint64_t rkn_retries_before = m_timing_rkn_rejected_trials;
  const std::uint64_t rkn_failures_before = m_timing_rkn_failures;
  const int run_number = get_run_number(topNode);
  const int event_number = get_event_number(topNode);
  Vec3 primary_vertex = m_fixed_primary_vertex;
  Tpc_PolyTrackVertexContainer *pattern_vertices = nullptr;
  std::map<int, Tracklet> tracklet_map;

  const auto track_build_start = std::chrono::steady_clock::now();
  if (m_use_pattern_cluster_tracks)
  {
    auto *clusters = findNode::getClass<Tpc_PolyClusterContainer>(
        topNode, m_tpc_sa_cluster_node);
    auto *tracks = findNode::getClass<Tpc_PolyTrackContainer>(
        topNode, m_tpc_sa_track_node);
    pattern_vertices = findNode::getClass<Tpc_PolyTrackVertexContainer>(
        topNode, m_tpc_sa_track_vertex_node);

    if (!clusters || !tracks)
    {
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << Name() << ": missing pattern-reco nodes "
                  << m_tpc_sa_cluster_node << "/" << m_tpc_sa_track_node << std::endl;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }

    if (!pattern_vertices && Verbosity() > 1)
    {
      std::cout << PHWHERE << Name() << ": missing pattern-reco vertex node "
                << m_tpc_sa_track_vertex_node
                << "; trackTree vertex_z will use the configured fallback vertex" << std::endl;
    }
    tracklet_map = build_pattern_tracklets(clusters, tracks);
  }
  else
  {
    auto *truth_points = findNode::getClass<PHG4HitContainer>(topNode, m_truth_point_node);
    if (!truth_points)
    {
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << Name() << ": missing truth point node "
                  << m_truth_point_node << std::endl;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }

    auto *truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, m_truth_info_node);
    if (!truth_info && Verbosity() > 0)
    {
      std::cout << PHWHERE << Name() << ": missing truth info node "
                << m_truth_info_node << "; no charged truth tracklets can be built" << std::endl;
    }

    primary_vertex = get_primary_vertex(truth_info);
    tracklet_map = build_tracklets(truth_points, truth_info);
  }
  const double track_build_seconds = std::chrono::duration<double>(
                                         std::chrono::steady_clock::now() - track_build_start)
                                         .count();
  if (m_print_timing)
  {
    std::cout << "[V0TimingStage] run=" << run_number
              << " event=" << event_number
              << " stage=track_build_done"
              << " tracks=" << tracklet_map.size()
              << " stage_s=" << track_build_seconds
              << " kalman_fit_s=" << (m_timing_kalman_fit_seconds - kalman_fit_before)
              << " rkn_propagations=" << (m_timing_rkn_propagations - rkn_propagations_before)
              << " rkn_steps=" << (m_timing_rkn_accepted_steps - rkn_steps_before)
              << std::endl;
  }

  const auto dca_cache_start = std::chrono::steady_clock::now();
  for (auto &entry : tracklet_map)
  {
    auto &tracklet = entry.second;
    tracklet.has_beamline_pca = track_pca_to_xy(
        tracklet, Vec3{0.0, 0.0, 0.0}, tracklet.beamline_pca, tracklet.rdca_zero);
    tracklet.has_pattern_vertex = choose_pattern_collision_vertex(
        tracklet, pattern_vertices, tracklet.pattern_vertex,
        tracklet.pattern_vertex_z_rms, tracklet.pattern_vertex_ntracks);
    const Vec3 &dca_vertex = tracklet.has_pattern_vertex
                                 ? tracklet.pattern_vertex
                                 : primary_vertex;
    tracklet.vertex_dca = fitted_track_dca_to_vertex(tracklet, dca_vertex);
    tracklet.has_vertex_dca =
        std::isfinite(tracklet.vertex_dca.first) &&
        std::isfinite(tracklet.vertex_dca.second);
  }
  const double dca_cache_seconds = std::chrono::duration<double>(
                                       std::chrono::steady_clock::now() - dca_cache_start)
                                       .count();
  m_timing_dca_cache_seconds += dca_cache_seconds;
  if (m_print_timing)
  {
    std::cout << "[V0TimingStage] run=" << run_number
              << " event=" << event_number
              << " stage=dca_cache_done"
              << " tracks=" << tracklet_map.size()
              << " stage_s=" << dca_cache_seconds
              << std::endl;
  }

  std::vector<const Tracklet *> tracklets;
  tracklets.reserve(tracklet_map.size());
  for (const auto &entry : tracklet_map)
  {
    tracklets.push_back(&entry.second);
  }

  const auto track_qa_start = std::chrono::steady_clock::now();
  for (const auto *tracklet : tracklets)
  {
    fill_track_row(*tracklet, primary_vertex, run_number, event_number);
    if (m_cluster_residual_tree)
    {
      fill_cluster_residual_rows(*tracklet, primary_vertex, run_number, event_number);
    }
  }
  const double track_qa_seconds = std::chrono::duration<double>(
                                      std::chrono::steady_clock::now() - track_qa_start)
                                      .count();
  if (m_print_timing)
  {
    std::cout << "[V0TimingStage] run=" << run_number
              << " event=" << event_number
              << " stage=track_qa_done"
              << " stage_s=" << track_qa_seconds
              << std::endl;
  }

  const auto pair_loop_start = std::chrono::steady_clock::now();
  std::uint64_t event_pairs_processed = 0;
  const std::uint64_t event_pairs_total =
      tracklets.size() > 1
          ? static_cast<std::uint64_t>(tracklets.size()) *
                static_cast<std::uint64_t>(tracklets.size() - 1) / 2
          : 0;
  if (m_print_timing)
  {
    std::cout << "[V0TimingStage] run=" << run_number
              << " event=" << event_number
              << " stage=pair_loop_start"
              << " pairs=" << event_pairs_total
              << std::endl;
  }
  for (std::size_t i = 0; i < tracklets.size(); ++i)
  {
    for (std::size_t j = i + 1; j < tracklets.size(); ++j)
    {
      make_pair_row(*tracklets[i], *tracklets[j], primary_vertex, run_number, event_number);
      ++event_pairs_processed;
      if (m_print_timing &&
          (event_pairs_processed == 1 || event_pairs_processed % 1000 == 0))
      {
        std::cout << "[V0TimingPair] run=" << run_number
                  << " event=" << event_number
                  << " done=" << event_pairs_processed
                  << " total=" << event_pairs_total
                  << " pair_loop_s=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - pair_loop_start).count()
                  << " kalman_pca_s=" << (m_timing_kalman_pca_seconds - kalman_pca_before)
                  << std::endl;
      }
    }
  }
  const double pair_loop_seconds = std::chrono::duration<double>(
                                       std::chrono::steady_clock::now() - pair_loop_start)
                                       .count();
  const double total_seconds = std::chrono::duration<double>(
                                   std::chrono::steady_clock::now() - event_start)
                                   .count();
  const double kalman_fit_seconds = m_timing_kalman_fit_seconds - kalman_fit_before;
  const double kalman_pca_seconds = m_timing_kalman_pca_seconds - kalman_pca_before;

  ++m_timing_events;
  m_timing_total_seconds += total_seconds;
  m_timing_track_build_seconds += track_build_seconds;
  m_timing_track_qa_seconds += track_qa_seconds;
  m_timing_pair_loop_seconds += pair_loop_seconds;

  if (m_print_timing)
  {
    std::cout << "[V0Timing] run=" << run_number
              << " event=" << event_number
              << " tracks=" << tracklets.size()
              << " total_s=" << total_seconds
              << " track_build_s=" << track_build_seconds
              << " kalman_fit_s=" << kalman_fit_seconds
              << " track_qa_s=" << track_qa_seconds
              << " pair_loop_s=" << pair_loop_seconds
              << " kalman_pca_s=" << kalman_pca_seconds
              << " rkn_propagations=" << (m_timing_rkn_propagations - rkn_propagations_before)
              << " rkn_steps=" << (m_timing_rkn_accepted_steps - rkn_steps_before)
              << " rkn_retries=" << (m_timing_rkn_rejected_trials - rkn_retries_before)
              << " rkn_failures=" << (m_timing_rkn_failures - rkn_failures_before)
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcV0CandidateTree::End(PHCompositeNode * /*topNode*/)
{
  if (m_file)
  {
    m_file->cd();
    if (m_pair_tree)
    {
      m_pair_tree->Write();
    }
    if (m_track_tree)
    {
      m_track_tree->Write();
    }
    if (m_cluster_residual_tree)
    {
      m_cluster_residual_tree->Write();
    }
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << ": pair counters: raw=" << m_counter_raw_pairs
              << " reject_charge=" << m_counter_reject_charge
              << " reject_preselection=" << m_counter_reject_preselection
              << " reject_pca=" << m_counter_reject_pca
              << " reject_pointing=" << m_counter_reject_pointing
              << " reject_ap=" << m_counter_reject_ap
              << " reject_pair_selection=" << m_counter_reject_pair_selection
              << " written=" << m_counter_written
              << " tracks_written=" << m_counter_tracks_written
              << " reject_helix_anchor=" << m_counter_reject_helix_anchor
              << " cluster_residuals_written=" << m_counter_cluster_residuals_written
              << std::endl;
  }

  if (m_print_timing)
  {
    std::cout << "[V0TimingSummary] events=" << m_timing_events
              << " total_s=" << m_timing_total_seconds
              << " track_build_s=" << m_timing_track_build_seconds
              << " kalman_fits=" << m_timing_kalman_fits
              << " kalman_fit_s=" << m_timing_kalman_fit_seconds
              << " rkn_s=" << m_timing_rkn_seconds
              << " dca_cache_s=" << m_timing_dca_cache_seconds
              << " track_qa_s=" << m_timing_track_qa_seconds
              << " pair_loop_s=" << m_timing_pair_loop_seconds
              << " kalman_pca_s=" << m_timing_kalman_pca_seconds
              << " rkn_propagations=" << m_timing_rkn_propagations
              << " rkn_steps=" << m_timing_rkn_accepted_steps
              << " rkn_retries=" << m_timing_rkn_rejected_trials
              << " rkn_failures=" << m_timing_rkn_failures
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcV0CandidateTree::get_event_number(PHCompositeNode *topNode) const
{
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    return event_header->get_EvtSequence();
  }
  if (auto *g4_event_header = findNode::getClass<PHG4EventHeader>(topNode, "EventHeader"))
  {
    return g4_event_header->get_EvtSequence();
  }
  return 0;
}

int TpcV0CandidateTree::get_run_number(PHCompositeNode *topNode) const
{
  if (auto *event_header = findNode::getClass<EventHeader>(topNode, "EventHeader"))
  {
    return event_header->get_RunNumber();
  }
  return 1;
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::get_primary_vertex(PHG4TruthInfoContainer *truth_info) const
{
  if (!m_use_truth_primary_vertex || !truth_info)
  {
    return m_fixed_primary_vertex;
  }

  const auto vtx_range = truth_info->GetPrimaryVtxRange();
  for (auto iter = vtx_range.first; iter != vtx_range.second; ++iter)
  {
    const auto *vtx = iter->second;
    if (vtx)
    {
      return {vtx->get_x(), vtx->get_y(), vtx->get_z()};
    }
  }

  return m_fixed_primary_vertex;
}

std::map<int, TpcV0CandidateTree::Tracklet> TpcV0CandidateTree::build_tracklets(
    PHG4HitContainer *truth_points,
    PHG4TruthInfoContainer *truth_info) const
{
  std::map<int, Tracklet> tracklets;
  if (!truth_points)
  {
    return tracklets;
  }

  const auto hit_range = truth_points->getHits();
  for (auto hit_iter = hit_range.first; hit_iter != hit_range.second; ++hit_iter)
  {
    const PHG4Hit *hit = hit_iter->second;
    if (!hit)
    {
      continue;
    }

    const int track_id = hit->get_trkid();
    Tracklet &tracklet = tracklets[track_id];
    if (tracklet.points.empty())
    {
      tracklet.track_id = track_id;
      tracklet.shower_id = hit->get_shower_id();
    }

    TruthPoint point;
    point.track_id = track_id;
    point.shower_id = hit->get_shower_id();
    point.layer = static_cast<int>(hit->get_layer());
    point.position = {hit->get_x(0), hit->get_y(0), hit->get_z(0)};
    point.momentum = {hit->get_px(0), hit->get_py(0), hit->get_pz(0)};
    point.t = hit->get_t(0);
    point.path = hit->get_path_length();
    if (finite(point.position) && finite(point.momentum))
    {
      tracklet.points.push_back(point);
    }
  }

  for (auto iter = tracklets.begin(); iter != tracklets.end();)
  {
    Tracklet &tracklet = iter->second;
    order_track_points(tracklet.points, m_point_order);
    tracklet.npoints = static_cast<int>(tracklet.points.size());
    tracklet.ntpc_clusters = static_cast<unsigned int>(tracklet.points.size());
    if (!tracklet.points.empty())
    {
      tracklet.side = tracklet.points.front().position.z < 0.0 ? 0 : 1;
    }

    if (tracklet.npoints < m_min_points)
    {
      iter = tracklets.erase(iter);
      continue;
    }

    if (truth_info)
    {
      const PHG4Particle *particle = truth_info->GetParticle(tracklet.track_id);
      if (particle)
      {
        tracklet.pid = particle->get_pid();
        tracklet.parent_id = particle->get_parent_id();
        tracklet.primary_id = particle->get_primary_id();
        tracklet.vtx_id = particle->get_vtx_id();
        tracklet.barcode = particle->get_barcode();
        tracklet.embed_id = truth_info->isEmbeded(tracklet.track_id);
        tracklet.is_primary = truth_info->is_primary(particle) ? 1 : 0;
        tracklet.charge = pdg_charge(tracklet.pid);
        tracklet.truth_momentum = {particle->get_px(), particle->get_py(), particle->get_pz()};
        tracklet.truth_e = particle->get_e();

        if (const PHG4Particle *parent = truth_info->GetParticle(tracklet.parent_id))
        {
          tracklet.parent_pid = parent->get_pid();
        }

        if (auto *vtx = truth_info->GetVtx(tracklet.vtx_id))
        {
          tracklet.truth_vertex = {vtx->get_x(), vtx->get_y(), vtx->get_z()};
          tracklet.truth_vt = vtx->get_t();
        }
      }
    }

    if (tracklet.charge == 0)
    {
      iter = tracklets.erase(iter);
      continue;
    }

    if (m_fit_kalman_tracks)
    {
      tracklet.has_kalman = fit_kalman(tracklet.points, tracklet.charge, tracklet.kalman);
      if (!tracklet.has_kalman || tracklet.kalman.states_smoothed.empty())
      {
        iter = tracklets.erase(iter);
        continue;
      }

      const auto state = TpcTrackKalmanFitter::propagation_state(tracklet.kalman, Vec3{});
      tracklet.position = TpcTrackKalmanFitter::state_position(state);
      tracklet.momentum = TpcTrackKalmanFitter::state_momentum(state);
    }
    else if (m_fit_helix_tracks)
    {
      tracklet.has_helix = fit_helix(tracklet.points, m_fit_first_points,
                                     tracklet.charge, m_bfield_t, tracklet.helix);
      if (!tracklet.has_helix)
      {
        iter = tracklets.erase(iter);
        continue;
      }
      tracklet.position = helix_point(tracklet.helix, tracklet.helix.theta_first);
      tracklet.momentum = helix_momentum(tracklet.helix, tracklet.helix.theta_first);
    }
    else
    {
      tracklet.position = tracklet.points.front().position;
      tracklet.momentum = tracklet.points.front().momentum;
    }

    assign_fit_quality(tracklet);

    if (!finite(tracklet.truth_momentum) || norm(tracklet.truth_momentum) <= 0.0)
    {
      tracklet.truth_momentum = tracklet.momentum;
    }

    ++iter;
  }

  return tracklets;
}

std::map<int, TpcV0CandidateTree::Tracklet> TpcV0CandidateTree::build_pattern_tracklets(
    Tpc_PolyClusterContainer *clusters,
    Tpc_PolyTrackContainer *tracks) const
{
  std::map<int, Tracklet> tracklets;

  if (!clusters || !tracks)
  {
    return tracklets;
  }

  std::map<unsigned int, std::vector<const Tpc_PolyCluster *>> clusters_by_source_id;
  for (unsigned int icluster = 0; icluster < clusters->size(); ++icluster)
  {
    const Tpc_PolyCluster *cluster = clusters->get_cluster(icluster);
    if (!cluster || !cluster->isValid())
    {
      continue;
    }
    clusters_by_source_id[cluster->get_source_assembled_track_id()].push_back(cluster);
  }

  for (unsigned int itrack = 0; itrack < tracks->size(); ++itrack)
  {
    const Tpc_PolyTrack *track = tracks->get_track(itrack);
    if (!track || !track->isValid() || track->get_fit_status() == 0)
    {
      continue;
    }

    const unsigned int source_id = track->get_source_assembled_track_id();
    const auto cluster_iter = clusters_by_source_id.find(source_id);
    if (cluster_iter == clusters_by_source_id.end())
    {
      continue;
    }

    const unsigned int pattern_track_id = track->get_track_id();
    const int track_id = static_cast<int>(pattern_track_id != 0 ? pattern_track_id : itrack + 1);

    Tracklet tracklet;
    tracklet.track_id = track_id;
    tracklet.shower_id = static_cast<int>(source_id);
    tracklet.charge = sign_to_charge(track->get_charge());
    tracklet.position = {track->get_x(), track->get_y(), track->get_z()};
    tracklet.momentum = {track->get_px(), track->get_py(), track->get_pz()};
    tracklet.truth_momentum = tracklet.momentum;
    tracklet.dedx = track->get_dedx();
    tracklet.has_dedx = std::isfinite(tracklet.dedx);

    std::map<int, unsigned int> side_counts;
    const auto &track_clusters = cluster_iter->second;
    tracklet.points.reserve(track_clusters.size());
    for (unsigned int icluster = 0; icluster < track_clusters.size(); ++icluster)
    {
      const Tpc_PolyCluster *cluster = track_clusters[icluster];
      if (!cluster || !cluster->isValid())
      {
        continue;
      }

      TruthPoint point;
      point.track_id = track_id;
      point.shower_id = static_cast<int>(source_id);
      point.layer = -1;
      if (cluster->size_hits() > 0)
      {
        point.layer = static_cast<int>(TrkrDefs::getLayer(cluster->get_hit_index(0).first));
      }
      point.position = {cluster->get_centroid_x(),
                        cluster->get_centroid_y(),
                        cluster->get_centroid_z()};
      point.momentum = tracklet.momentum;
      point.t = 0.0;
      point.path = point.layer >= 0 ? static_cast<double>(point.layer)
                                    : static_cast<double>(icluster);
      if (finite(point.position))
      {
        tracklet.points.push_back(point);
        ++side_counts[cluster->get_side()];
      }
    }

    if (!side_counts.empty())
    {
      tracklet.side = std::max_element(
                          side_counts.begin(), side_counts.end(),
                          [](const auto &lhs, const auto &rhs)
                          { return lhs.second < rhs.second; })
                          ->first;
    }
    tracklet.ntpc_clusters = track->get_nclusters() > 0
                                 ? track->get_nclusters()
                                 : static_cast<unsigned int>(tracklet.points.size());

    const bool has_upstream_state = finite(tracklet.position) &&
                                    finite(tracklet.momentum) &&
                                    norm(tracklet.momentum) > 0.0;
    if (!finalize_pattern_tracklet(tracklet, has_upstream_state))
    {
      continue;
    }

    tracklet.truth_momentum = tracklet.momentum;
    tracklets[track_id] = std::move(tracklet);
  }
  return tracklets;
}

bool TpcV0CandidateTree::finalize_pattern_tracklet(Tracklet &tracklet,
                                                   const bool has_upstream_state) const
{
  order_track_points(tracklet.points, m_point_order);
  tracklet.npoints = static_cast<int>(tracklet.points.size());
  if (tracklet.npoints < m_min_points || tracklet.charge == 0)
  {
    return false;
  }

  if (m_fit_kalman_tracks)
  {
    tracklet.has_kalman = fit_kalman(tracklet.points, tracklet.charge, tracklet.kalman);
    if (!tracklet.has_kalman || tracklet.kalman.states_smoothed.empty())
    {
      return false;
    }

    const auto state = TpcTrackKalmanFitter::propagation_state(tracklet.kalman, Vec3{});
    tracklet.position = TpcTrackKalmanFitter::state_position(state);
    tracklet.momentum = TpcTrackKalmanFitter::state_momentum(state);
  }
  else if (m_use_final_track_helix && has_upstream_state)
  {
    tracklet.has_helix = helix_from_state(tracklet.position, tracklet.momentum,
                                          tracklet.charge, m_bfield_t, tracklet.helix);
    if (!tracklet.has_helix)
    {
      return false;
    }

    tracklet.has_helix_search_range =
        TpcTrackHelixFitter::measurement_anchored_search_range(
            tracklet.helix, tracklet.points,
            m_final_track_helix_max_upstream_cm,
            m_final_track_helix_downstream_margin_cm,
            tracklet.helix_search_range);
    if (!tracklet.has_helix_search_range)
    {
      ++m_counter_reject_helix_anchor;
      return false;
    }

    tracklet.position = helix_point(tracklet.helix, tracklet.helix.theta_first);
    tracklet.momentum = helix_momentum(tracklet.helix, tracklet.helix.theta_first);
  }
  else if (m_fit_helix_tracks)
  {
    tracklet.has_helix = fit_helix(tracklet.points, m_fit_first_points,
                                   tracklet.charge, m_bfield_t, tracklet.helix);
    if (!tracklet.has_helix)
    {
      return false;
    }

    tracklet.position = helix_point(tracklet.helix, tracklet.helix.theta_first);
    tracklet.momentum = helix_momentum(tracklet.helix, tracklet.helix.theta_first);
  }
  else if (!has_upstream_state)
  {
    if (tracklet.points.size() < 2)
    {
      return false;
    }
    tracklet.position = tracklet.points.front().position;
    tracklet.momentum = subtract(tracklet.points[1].position,
                                 tracklet.points.front().position);
  }

  assign_fit_quality(tracklet);
  return finite(tracklet.position) && finite(tracklet.momentum) &&
         norm(tracklet.momentum) > 0.0;
}

bool TpcV0CandidateTree::track_pca_to_xy(const Tracklet &tracklet,
                                         const Vec3 &beamline,
                                         Vec3 &pca,
                                         double &signed_dca_xy) const
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  pca = {nan, nan, nan};
  signed_dca_xy = nan;

  if (tracklet.has_kalman && !tracklet.kalman.states_smoothed.empty())
  {
    const auto state = TpcTrackKalmanFitter::propagation_state(tracklet.kalman, beamline);
    const Vec3 position = TpcTrackKalmanFitter::state_position(state);
    const Vec3 momentum = TpcTrackKalmanFitter::state_momentum(state);
    const double qop_t = state[TpcTrackKalmanFitter::QOverPt];
    const double omega = 0.003 * tracklet.kalman.bfield_t * qop_t;
    if (std::abs(omega) >= 1.0e-10)
    {
      const double radius = std::abs(1.0 / omega);
      const double center_x = state[TpcTrackKalmanFitter::X] -
                              std::sin(state[TpcTrackKalmanFitter::Phi]) / omega;
      const double center_y = state[TpcTrackKalmanFitter::Y] +
                              std::cos(state[TpcTrackKalmanFitter::Phi]) / omega;
      const double dx = beamline.x - center_x;
      const double dy = beamline.y - center_y;
      const double center_distance = std::hypot(dx, dy);
      if (center_distance > 1.0e-12)
      {
        const double closest_x = center_x + radius * dx / center_distance;
        const double closest_y = center_y + radius * dy / center_distance;
        const double theta0 = std::atan2(position.y - center_y, position.x - center_x);
        const double theta_closest = std::atan2(closest_y - center_y,
                                                closest_x - center_x);
        const double path_cm = normalize_phi(theta_closest - theta0) / omega;
        TpcKalmanConfig config = m_kalman_config;
        config.bfield_t = tracklet.kalman.bfield_t;
        config.magnetic_field = tracklet.kalman.magnetic_field;
        config.analytic_uniform_propagation = tracklet.kalman.analytic_uniform_propagation;
        const auto closest_state = TpcTrackKalmanFitter::propagate_state(
            state, path_cm, config, tracklet.kalman.mass_gev);
        pca = TpcTrackKalmanFitter::state_position(closest_state);
        signed_dca_xy = center_distance - radius;
        return finite(pca) && std::isfinite(signed_dca_xy);
      }

      pca = position;
      signed_dca_xy = -radius;
      return finite(pca);
    }

    const double transverse_momentum2 = square(momentum.x) + square(momentum.y);
    if (transverse_momentum2 <= 0.0)
    {
      return false;
    }
    const Vec3 relative = subtract(position, beamline);
    const double scale_to_pca = -(relative.x * momentum.x + relative.y * momentum.y) /
                                transverse_momentum2;
    pca = add(position, scale(momentum, scale_to_pca));
    signed_dca_xy = (relative.x * momentum.y - relative.y * momentum.x) /
                    std::sqrt(transverse_momentum2);
    return finite(pca) && std::isfinite(signed_dca_xy);
  }

  if (tracklet.has_helix && tracklet.helix.radius > 0.0)
  {
    const double dx = beamline.x - tracklet.helix.cx;
    const double dy = beamline.y - tracklet.helix.cy;
    const double center_distance = std::hypot(dx, dy);
    double theta = tracklet.helix.theta_first;
    if (center_distance > 1.0e-12)
    {
      const double theta_raw = std::atan2(dy, dx);
      theta = unwrap_to_near(theta_raw, tracklet.helix.theta_first);
    }
    pca = helix_point(tracklet.helix, theta);
    signed_dca_xy = center_distance - tracklet.helix.radius;
    return finite(pca) && std::isfinite(signed_dca_xy);
  }

  const double transverse_momentum2 = square(tracklet.momentum.x) + square(tracklet.momentum.y);
  if (transverse_momentum2 <= 0.0)
  {
    return false;
  }
  const Vec3 relative = subtract(tracklet.position, beamline);
  const double scale_to_pca = -(relative.x * tracklet.momentum.x +
                                relative.y * tracklet.momentum.y) /
                              transverse_momentum2;
  pca = add(tracklet.position, scale(tracklet.momentum, scale_to_pca));
  signed_dca_xy = (relative.x * tracklet.momentum.y -
                   relative.y * tracklet.momentum.x) /
                  std::sqrt(transverse_momentum2);
  return finite(pca) && std::isfinite(signed_dca_xy);
}

bool TpcV0CandidateTree::choose_pattern_collision_vertex(
    const Tracklet &tracklet,
    Tpc_PolyTrackVertexContainer *vertices,
    Vec3 &vertex,
    double &z_rms,
    unsigned int &ntracks) const
{
  if (!vertices || vertices->get_collision_vertex_valid() == 0)
  {
    return false;
  }

  double best_dz = std::numeric_limits<double>::max();
  const unsigned int count = vertices->get_collision_vertex_count();
  for (unsigned int index = 0; index < count; ++index)
  {
    const Vec3 candidate{vertices->get_collision_x(index),
                         vertices->get_collision_y(index),
                         vertices->get_collision_z(index)};
    if (!finite(candidate))
    {
      continue;
    }

    Vec3 candidate_pca;
    double candidate_dca = 0.0;
    if (!track_pca_to_xy(tracklet, candidate, candidate_pca, candidate_dca))
    {
      continue;
    }
    const double dz = std::abs(candidate_pca.z - candidate.z);
    if (dz < best_dz)
    {
      best_dz = dz;
      vertex = candidate;
      z_rms = vertices->get_collision_z_rms(index);
      ntracks = vertices->get_collision_ntracks(index);
    }
  }
  return best_dz < std::numeric_limits<double>::max();
}

bool TpcV0CandidateTree::make_pair_row(const Tracklet &track1, const Tracklet &track2,
                                       const Vec3 &primary_vertex,
                                       const int run_number,
                                       const int event_number)
{
  ++m_counter_raw_pairs;

  if (!m_write_same_sign_pairs && track1.charge == track2.charge)
  {
    ++m_counter_reject_charge;
    return false;
  }

  if (!passes_preselection(track1, track2, primary_vertex))
  {
    ++m_counter_reject_preselection;
    return false;
  }

  Vec3 pca1;
  Vec3 pca2;
  Vec3 mom1 = track1.momentum;
  Vec3 mom2 = track2.momentum;
  double pair_dca = 0.0;
  double theta1 = quiet_nan();
  double theta2 = quiet_nan();
  std::pair<double, double> dca1;
  std::pair<double, double> dca2;

  if (m_fit_kalman_tracks && track1.has_kalman && track2.has_kalman)
  {
    const auto pca_start = std::chrono::steady_clock::now();
    auto candidates = kalman_pca_candidates(
        track1.kalman, track2.kalman, m_kalman_config, primary_vertex,
        m_kalman_max_upstream_cm, m_kalman_downstream_margin_cm,
        m_coarse_steps, m_pca_candidates);
    m_timing_kalman_pca_seconds += std::chrono::duration<double>(
                                       std::chrono::steady_clock::now() - pca_start)
                                       .count();
    if (candidates.empty())
    {
      ++m_counter_reject_pca;
      return false;
    }

    auto best = candidates.front();
    if (m_prefer_positive_pointing)
    {
      double best_score = std::numeric_limits<double>::max();
      for (const auto &candidate : candidates)
      {
        const Vec3 cand_mom1 = kalman_momentum(track1.kalman, candidate.s1, m_kalman_config, primary_vertex);
        const Vec3 cand_mom2 = kalman_momentum(track2.kalman, candidate.s2, m_kalman_config, primary_vertex);
        const Vec3 cand_vertex = scale(add(candidate.pca1, candidate.pca2), 0.5);
        const Vec3 flight = subtract(cand_vertex, primary_vertex);
        const Vec3 total_mom = add(cand_mom1, cand_mom2);
        const double cos_theta = vector_cosine(flight, total_mom);
        const double penalty = (std::isfinite(cos_theta) && cos_theta > 0.0) ? 0.0 : 1000.0;
        const double score = penalty + candidate.dca - 1e-3 * cos_theta;
        if (score < best_score)
        {
          best_score = score;
          best = candidate;
        }
      }
    }

    pca1 = best.pca1;
    pca2 = best.pca2;
    pair_dca = best.dca;
    theta1 = best.s1;
    theta2 = best.s2;
    mom1 = kalman_momentum(track1.kalman, theta1, m_kalman_config, primary_vertex);
    mom2 = kalman_momentum(track2.kalman, theta2, m_kalman_config, primary_vertex);
    dca1 = track1.has_vertex_dca
               ? track1.vertex_dca
               : TpcTrackKalmanFitter::dca_to_vertex(
                     track1.kalman, primary_vertex, &m_kalman_config);
    dca2 = track2.has_vertex_dca
               ? track2.vertex_dca
               : TpcTrackKalmanFitter::dca_to_vertex(
                     track2.kalman, primary_vertex, &m_kalman_config);
  }
  else if (m_fit_helix_tracks && track1.has_helix && track2.has_helix)
  {
    std::vector<HelixPca> candidates;
    if (track1.has_helix_search_range && track2.has_helix_search_range)
    {
      candidates = TpcTrackHelixFitter::pca_candidates_in_ranges(
          track1.helix, track2.helix,
          track1.helix_search_range, track2.helix_search_range,
          m_coarse_steps, m_pca_candidates);
    }
    else
    {
      candidates = helix_helix_pca_candidates(
          track1.helix, track2.helix, m_theta_extension, m_coarse_steps,
          m_downstream_margin, m_pca_candidates);
    }
    if (candidates.empty())
    {
      ++m_counter_reject_pca;
      return false;
    }

    auto best = candidates.front();
    if (m_prefer_positive_pointing)
    {
      double best_score = std::numeric_limits<double>::max();
      for (const auto &candidate : candidates)
      {
        const Vec3 cand_mom1 = helix_momentum(track1.helix, candidate.theta1);
        const Vec3 cand_mom2 = helix_momentum(track2.helix, candidate.theta2);
        const Vec3 cand_vertex = scale(add(candidate.pca1, candidate.pca2), 0.5);
        const Vec3 flight = subtract(cand_vertex, primary_vertex);
        const Vec3 total_mom = add(cand_mom1, cand_mom2);
        const double cos_theta = vector_cosine(flight, total_mom);
        const double penalty = (std::isfinite(cos_theta) && cos_theta > 0.0) ? 0.0 : 1000.0;
        const double score = penalty + candidate.dca - 1e-3 * cos_theta;
        if (score < best_score)
        {
          best_score = score;
          best = candidate;
        }
      }
    }

    pca1 = best.pca1;
    pca2 = best.pca2;
    pair_dca = best.dca;
    theta1 = best.theta1;
    theta2 = best.theta2;
    mom1 = helix_momentum(track1.helix, theta1);
    mom2 = helix_momentum(track2.helix, theta2);
    dca1 = helix_dca_to_vertex(track1.helix, primary_vertex);
    dca2 = helix_dca_to_vertex(track2.helix, primary_vertex);
  }
  else
  {
    LinePca pca;
    if (!line_line_pca(track1.position, track1.momentum, track2.position, track2.momentum, pca, true))
    {
      ++m_counter_reject_pca;
      return false;
    }
    pca1 = pca.pca1;
    pca2 = pca.pca2;
    pair_dca = pca.dca;
    dca1 = track_dca_to_vertex(track1.position, track1.momentum, primary_vertex);
    dca2 = track_dca_to_vertex(track2.position, track2.momentum, primary_vertex);
  }

  const Vec3 pair_vertex = scale(add(pca1, pca2), 0.5);
  const Vec3 total_mom = add(mom1, mom2);
  const Vec3 flight = subtract(pair_vertex, primary_vertex);
  const double cos_theta = vector_cosine(flight, total_mom);
  if (!std::isfinite(cos_theta) || norm(flight) <= 0.0 || norm(total_mom) <= 0.0)
  {
    ++m_counter_reject_pointing;
    return false;
  }

  const Vec3 &pplus = (track1.charge > 0) ? mom1 : mom2;
  const Vec3 &pminus = (track1.charge > 0) ? mom2 : mom1;
  double alpha = 0.0;
  double qt = 0.0;
  if (!armenteros(pplus, pminus, alpha, qt))
  {
    ++m_counter_reject_ap;
    return false;
  }

  if (!passes_pair_selection(pca1, pca2, pair_vertex, primary_vertex,
                             pair_dca, cos_theta, alpha))
  {
    ++m_counter_reject_pair_selection;
    return false;
  }

  const Vec3 &truth_pplus = (track1.charge > 0) ? track1.truth_momentum : track2.truth_momentum;
  const Vec3 &truth_pminus = (track1.charge > 0) ? track2.truth_momentum : track1.truth_momentum;
  double truth_alpha = quiet_nan();
  double truth_qt = quiet_nan();
  armenteros(truth_pplus, truth_pminus, truth_alpha, truth_qt);

  Vec3 true_decay{quiet_nan(), quiet_nan(), quiet_nan()};
  double pca_to_true_3d = quiet_nan();
  double pca_to_true_xy = quiet_nan();
  double pca_to_true_z = quiet_nan();
  if (track1.parent_id != 0 && track1.parent_id == track2.parent_id)
  {
    true_decay = scale(add(track1.truth_vertex, track2.truth_vertex), 0.5);
    const Vec3 delta = subtract(pair_vertex, true_decay);
    pca_to_true_3d = norm(delta);
    pca_to_true_xy = std::sqrt(square(delta.x) + square(delta.y));
    pca_to_true_z = std::abs(delta.z);
  }

  const Vec3 positive_mom = (track1.charge > 0) ? mom1 : mom2;
  const Vec3 negative_mom = (track1.charge > 0) ? mom2 : mom1;

  reset_pair_row();
  m_pair.run = run_number;
  m_pair.evt = event_number;
  m_pair.cross1 = 0;
  m_pair.cross2 = 0;
  m_pair.px1 = static_cast<float>(mom1.x);
  m_pair.py1 = static_cast<float>(mom1.y);
  m_pair.pz1 = static_cast<float>(mom1.z);
  m_pair.px2 = static_cast<float>(mom2.x);
  m_pair.py2 = static_cast<float>(mom2.y);
  m_pair.pz2 = static_cast<float>(mom2.z);
  m_pair.dca_xy1 = static_cast<float>(dca1.first);
  m_pair.dca_z1 = static_cast<float>(dca1.second);
  m_pair.dca_xy2 = static_cast<float>(dca2.first);
  m_pair.dca_z2 = static_cast<float>(dca2.second);
  m_pair.pairDCA = static_cast<float>(pair_dca);
  m_pair.alpha = static_cast<float>(alpha);
  m_pair.qT = static_cast<float>(qt);
  m_pair.charge1 = static_cast<float>(track1.charge);
  m_pair.charge2 = static_cast<float>(track2.charge);
  m_pair.dedx_1 = track1.has_dedx ? static_cast<float>(track1.dedx) : quiet_nan();
  m_pair.dedx_2 = track2.has_dedx ? static_cast<float>(track2.dedx) : quiet_nan();
  m_pair.cosThetaReco = static_cast<float>(cos_theta);
  m_pair.Lproj = static_cast<float>(norm(flight));

  m_pair.pca_x = static_cast<float>(pair_vertex.x);
  m_pair.pca_y = static_cast<float>(pair_vertex.y);
  m_pair.pca_z = static_cast<float>(pair_vertex.z);
  m_pair.pca1_x = static_cast<float>(pca1.x);
  m_pair.pca1_y = static_cast<float>(pca1.y);
  m_pair.pca1_z = static_cast<float>(pca1.z);
  m_pair.pca2_x = static_cast<float>(pca2.x);
  m_pair.pca2_y = static_cast<float>(pca2.y);
  m_pair.pca2_z = static_cast<float>(pca2.z);

  m_pair.v0_px = static_cast<float>(total_mom.x);
  m_pair.v0_py = static_cast<float>(total_mom.y);
  m_pair.v0_pz = static_cast<float>(total_mom.z);
  m_pair.v0_pt = static_cast<float>(pt(total_mom));
  m_pair.mass_Kshort = static_cast<float>(invariant_mass(mom1, kPionMass, mom2, kPionMass));
  m_pair.mass_Lambda = static_cast<float>(invariant_mass(positive_mom, kProtonMass, negative_mom, kPionMass));
  m_pair.mass_AntiLambda = static_cast<float>(invariant_mass(positive_mom, kPionMass, negative_mom, kProtonMass));

  m_pair.true_decay_x = static_cast<float>(true_decay.x);
  m_pair.true_decay_y = static_cast<float>(true_decay.y);
  m_pair.true_decay_z = static_cast<float>(true_decay.z);
  m_pair.pca_to_true_3d = static_cast<float>(pca_to_true_3d);
  m_pair.pca_to_true_xy = static_cast<float>(pca_to_true_xy);
  m_pair.pca_to_true_z = static_cast<float>(pca_to_true_z);
  m_pair.truth_alpha = static_cast<float>(truth_alpha);
  m_pair.truth_qT = static_cast<float>(truth_qt);
  m_pair.delta_alpha = static_cast<float>(alpha - truth_alpha);
  m_pair.delta_qT = static_cast<float>(qt - truth_qt);
  m_pair.truth_px1 = static_cast<float>(track1.truth_momentum.x);
  m_pair.truth_py1 = static_cast<float>(track1.truth_momentum.y);
  m_pair.truth_pz1 = static_cast<float>(track1.truth_momentum.z);
  m_pair.truth_px2 = static_cast<float>(track2.truth_momentum.x);
  m_pair.truth_py2 = static_cast<float>(track2.truth_momentum.y);
  m_pair.truth_pz2 = static_cast<float>(track2.truth_momentum.z);
  m_pair.cos_mom1_truth = static_cast<float>(vector_cosine(mom1, track1.truth_momentum));
  m_pair.cos_mom2_truth = static_cast<float>(vector_cosine(mom2, track2.truth_momentum));
  m_pair.pca_theta1 = static_cast<float>(theta1);
  m_pair.pca_theta2 = static_cast<float>(theta2);
  if (track1.has_kalman)
  {
    m_pair.kalman_chi2_1 = static_cast<float>(track1.kalman.chi2);
    m_pair.kalman_ndof1 = track1.kalman.ndof;
    m_pair.kalman_chi2_ndf1 =
        (track1.kalman.ndof > 0 && std::isfinite(track1.kalman.chi2))
            ? static_cast<float>(track1.kalman.chi2 / track1.kalman.ndof)
            : quiet_nan();
  }
  if (track2.has_kalman)
  {
    m_pair.kalman_chi2_2 = static_cast<float>(track2.kalman.chi2);
    m_pair.kalman_ndof2 = track2.kalman.ndof;
    m_pair.kalman_chi2_ndf2 =
        (track2.kalman.ndof > 0 && std::isfinite(track2.kalman.chi2))
            ? static_cast<float>(track2.kalman.chi2 / track2.kalman.ndof)
            : quiet_nan();
  }
  m_pair.quality1 = std::isfinite(track1.fit_chi2_ndf)
                        ? static_cast<float>(track1.fit_chi2_ndf)
                        : quiet_nan();
  m_pair.quality2 = std::isfinite(track2.fit_chi2_ndf)
                        ? static_cast<float>(track2.fit_chi2_ndf)
                        : quiet_nan();
  m_pair.track_id1 = track1.track_id;
  m_pair.track_id2 = track2.track_id;
  m_pair.pid1 = track1.pid;
  m_pair.pid2 = track2.pid;
  m_pair.parent_id1 = track1.parent_id;
  m_pair.parent_id2 = track2.parent_id;
  m_pair.parent_pid = (track1.parent_id != 0 && track1.parent_id == track2.parent_id) ? track1.parent_pid : 0;
  m_pair.npoints1 = static_cast<short>(track1.npoints);
  m_pair.npoints2 = static_cast<short>(track2.npoints);

  m_pair_tree->Fill();
  ++m_counter_written;
  return true;
}

void TpcV0CandidateTree::fill_track_row(const Tracklet &tracklet,
                                        const Vec3 &primary_vertex,
                                        const int run_number,
                                        const int event_number)
{
  reset_track_row();
  const float nan = quiet_nan();

  m_track.run = run_number;
  m_track.evt = event_number;
  m_track.track_id = tracklet.track_id;
  m_track.shower_id = tracklet.shower_id;
  m_track.pid = tracklet.pid;
  m_track.parent_id = tracklet.parent_id;
  m_track.parent_pid = tracklet.parent_pid;
  m_track.charge = static_cast<double>(tracklet.charge);
  m_track.side = tracklet.side;
  m_track.npoints = tracklet.npoints;
  m_track.ntpc_clusters = tracklet.ntpc_clusters;
  m_track.has_helix = tracklet.has_helix ? 1 : 0;
  m_track.has_kalman = tracklet.has_kalman ? 1 : 0;
  m_track.is_primary = tracklet.is_primary;

  if (m_kalman_config.collect_innovation_components)
  {
    m_track.kalman_measurement_sigma_r =
        static_cast<float>(m_kalman_config.meas_sigma_r_cm);
    m_track.kalman_measurement_sigma_rphi =
        static_cast<float>(m_kalman_config.meas_sigma_rphi_cm);
    m_track.kalman_measurement_sigma_z =
        static_cast<float>(m_kalman_config.meas_sigma_z_cm);
  }

  Vec3 row_position = tracklet.position;
  Vec3 row_momentum = tracklet.momentum;
  std::array<double, TpcTrackKalmanFitter::StateDim> kalman_row_state{};
  bool has_kalman_row_state = false;
  if (tracklet.has_kalman)
  {
    kalman_row_state = TpcTrackKalmanFitter::propagation_state(tracklet.kalman, primary_vertex);
    row_position = TpcTrackKalmanFitter::state_position(kalman_row_state);
    row_momentum = TpcTrackKalmanFitter::state_momentum(kalman_row_state);
    has_kalman_row_state = true;
  }

  m_track.px = row_momentum.x;
  m_track.py = row_momentum.y;
  m_track.pz = row_momentum.z;
  m_track.pt = pt(row_momentum);
  m_track.p = norm(row_momentum);
  m_track.eta = m_track.pt > 0.0
                    ? std::asinh(row_momentum.z / m_track.pt)
                    : static_cast<double>(nan);
  m_track.dedx = tracklet.has_dedx ? tracklet.dedx : static_cast<double>(nan);
  m_track.x = static_cast<float>(row_position.x);
  m_track.y = static_cast<float>(row_position.y);
  m_track.z = static_cast<float>(row_position.z);

  if (!tracklet.points.empty())
  {
    const auto &first = tracklet.points.front().position;
    const auto &last = tracklet.points.back().position;
    m_track.first_x = static_cast<float>(first.x);
    m_track.first_y = static_cast<float>(first.y);
    m_track.first_z = static_cast<float>(first.z);
    m_track.first_r = static_cast<float>(pt(first));
    m_track.last_x = static_cast<float>(last.x);
    m_track.last_y = static_cast<float>(last.y);
    m_track.last_z = static_cast<float>(last.z);
    m_track.last_r = static_cast<float>(pt(last));
  }

  m_track.cluster_index.reserve(tracklet.points.size());
  m_track.cluster_side.reserve(tracklet.points.size());
  m_track.layer.reserve(tracklet.points.size());
  m_track.cluster_z.reserve(tracklet.points.size());
  m_track.cluster_r.reserve(tracklet.points.size());
  m_track.cluster_phi.reserve(tracklet.points.size());
  m_track.residual_z.reserve(tracklet.points.size());
  m_track.residual_r.reserve(tracklet.points.size());
  m_track.residual_rphi.reserve(tracklet.points.size());

  double previous_theta = tracklet.has_helix ? tracklet.helix.theta_first : 0.0;
  bool have_previous_theta = false;
  for (std::size_t index = 0; index < tracklet.points.size(); ++index)
  {
    const auto &point = tracklet.points[index];
    const Vec3 &cluster = point.position;
    const double cluster_r = pt(cluster);
    const double cluster_phi = std::atan2(cluster.y, cluster.x);

    double residual_r = std::numeric_limits<double>::quiet_NaN();
    double residual_rphi = std::numeric_limits<double>::quiet_NaN();
    double residual_z = std::numeric_limits<double>::quiet_NaN();

    if (tracklet.has_kalman && index < tracklet.kalman.states_smoothed.size())
    {
      const Vec3 fit_position = TpcTrackKalmanFitter::state_position(tracklet.kalman.states_smoothed[index]);
      const Vec3 delta = subtract(cluster, fit_position);
      const double fit_phi = std::atan2(fit_position.y, fit_position.x);
      residual_r = std::cos(fit_phi) * delta.x + std::sin(fit_phi) * delta.y;
      residual_rphi = -std::sin(fit_phi) * delta.x + std::cos(fit_phi) * delta.y;
      residual_z = delta.z;
    }
    else if (tracklet.has_helix)
    {
      double theta_hint = std::atan2(cluster.y - tracklet.helix.cy,
                                     cluster.x - tracklet.helix.cx);
      theta_hint = unwrap_to_near(theta_hint,
                                  have_previous_theta ? previous_theta
                                                      : tracklet.helix.theta_first);

      double theta = theta_hint;
      Vec3 fit_position;
      if (!helix_point_at_beam_radius(tracklet.helix, cluster_r, theta_hint, theta, fit_position))
      {
        fit_position = helix_point(tracklet.helix, theta_hint);
        theta = theta_hint;
      }

      if (finite(fit_position))
      {
        previous_theta = theta;
        have_previous_theta = true;
        const double fit_r = pt(fit_position);
        const double fit_phi = std::atan2(fit_position.y, fit_position.x);
        residual_r = cluster_r - fit_r;
        residual_rphi = cluster_r * normalize_phi(cluster_phi - fit_phi);
        residual_z = cluster.z - fit_position.z;
      }
    }

    m_track.cluster_index.push_back(static_cast<unsigned int>(index));
    m_track.cluster_side.push_back(tracklet.side);
    m_track.layer.push_back(static_cast<unsigned int>(std::max(point.layer, 0)));
    m_track.cluster_z.push_back(cluster.z);
    m_track.cluster_r.push_back(cluster_r);
    m_track.cluster_phi.push_back(cluster_phi);
    m_track.residual_z.push_back(std::isfinite(residual_z) ? residual_z : static_cast<double>(nan));
    m_track.residual_r.push_back(std::isfinite(residual_r) ? residual_r : static_cast<double>(nan));
    m_track.residual_rphi.push_back(std::isfinite(residual_rphi) ? residual_rphi : static_cast<double>(nan));
  }

  const Vec3 &track_vertex = tracklet.has_pattern_vertex
                                 ? tracklet.pattern_vertex
                                 : primary_vertex;
  auto dca = tracklet.vertex_dca;
  if (!tracklet.has_vertex_dca)
  {
    dca = fitted_track_dca_to_vertex(tracklet, track_vertex);
  }
  m_track.dca_xy = static_cast<float>(dca.first);
  m_track.dca_z = static_cast<float>(dca.second);
  m_track.vertex_x = track_vertex.x;
  m_track.vertex_y = track_vertex.y;
  m_track.vertex_z = track_vertex.z;
  m_track.vertex_from_upstream = tracklet.has_pattern_vertex ? 1 : 0;
  if (tracklet.has_pattern_vertex)
  {
    m_track.vertex_z_rms = tracklet.pattern_vertex_z_rms;
    m_track.vertex_ntracks = tracklet.pattern_vertex_ntracks;
  }
  if (tracklet.has_beamline_pca)
  {
    m_track.pca_x = tracklet.beamline_pca.x;
    m_track.pca_y = tracklet.beamline_pca.y;
    m_track.pca_z = tracklet.beamline_pca.z;
    m_track.rDCA_zero = tracklet.rdca_zero;
    m_track.zDCA = tracklet.beamline_pca.z - track_vertex.z;
  }

  if (tracklet.has_helix)
  {
    m_track.helix_cx = static_cast<float>(tracklet.helix.cx);
    m_track.helix_cy = static_cast<float>(tracklet.helix.cy);
    m_track.helix_radius = static_cast<float>(tracklet.helix.radius);
    m_track.helix_z0 = static_cast<float>(tracklet.helix.z0);
    m_track.helix_pitch = static_cast<float>(tracklet.helix.pitch);
    m_track.helix_theta_first = static_cast<float>(tracklet.helix.theta_first);
    m_track.helix_theta_last = static_cast<float>(tracklet.helix.theta_last);
    m_track.helix_direction = static_cast<float>(tracklet.helix.direction);
    if (tracklet.has_helix_search_range)
    {
      const auto &range = tracklet.helix_search_range;
      m_track.helix_search_anchored = 1;
      m_track.helix_anchor_point_index = range.anchor_point_index;
      m_track.helix_anchor_theta = static_cast<float>(range.anchor_theta);
      m_track.helix_anchor_path_cm = static_cast<float>(range.anchor_path_cm);
      m_track.helix_anchor_residual_cm = static_cast<float>(range.anchor_residual_cm);
      m_track.helix_search_theta_min = static_cast<float>(range.theta_min);
      m_track.helix_search_theta_max = static_cast<float>(range.theta_max);
      m_track.helix_search_upstream_cm = static_cast<float>(range.upstream_cm);
      m_track.helix_search_downstream_cm = static_cast<float>(range.downstream_cm);
    }
  }

  if (tracklet.has_kalman)
  {
    m_track.kalman_chi2 = static_cast<float>(tracklet.kalman.chi2);
    m_track.kalman_ndof = tracklet.kalman.ndof;
    m_track.kalman_naccepted = static_cast<unsigned int>(tracklet.kalman.naccepted);
    m_track.kalman_nrejected = static_cast<unsigned int>(tracklet.kalman.nrejected);
    m_track.kalman_measurement_chi2 = tracklet.kalman.measurement_chi2;
    m_track.kalman_measurement_used = tracklet.kalman.measurement_used;
    if (m_kalman_config.collect_innovation_components)
    {
      m_track.kalman_measurement_in_seed = tracklet.kalman.measurement_in_seed;
      m_track.kalman_innovation_residual_r = tracklet.kalman.innovation_residual_r;
      m_track.kalman_innovation_residual_rphi = tracklet.kalman.innovation_residual_rphi;
      m_track.kalman_innovation_residual_z = tracklet.kalman.innovation_residual_z;
      m_track.kalman_prediction_sigma_r = tracklet.kalman.prediction_sigma_r;
      m_track.kalman_prediction_sigma_rphi = tracklet.kalman.prediction_sigma_rphi;
      m_track.kalman_prediction_sigma_z = tracklet.kalman.prediction_sigma_z;
      m_track.kalman_innovation_sigma_r = tracklet.kalman.innovation_sigma_r;
      m_track.kalman_innovation_sigma_rphi = tracklet.kalman.innovation_sigma_rphi;
      m_track.kalman_innovation_sigma_z = tracklet.kalman.innovation_sigma_z;
      m_track.kalman_innovation_rho_r_rphi = tracklet.kalman.innovation_rho_r_rphi;
      m_track.kalman_innovation_rho_r_z = tracklet.kalman.innovation_rho_r_z;
      m_track.kalman_innovation_rho_rphi_z = tracklet.kalman.innovation_rho_rphi_z;
      m_track.kalman_innovation_whitened_0 = tracklet.kalman.innovation_whitened_0;
      m_track.kalman_innovation_whitened_1 = tracklet.kalman.innovation_whitened_1;
      m_track.kalman_innovation_whitened_2 = tracklet.kalman.innovation_whitened_2;
    }
    if (has_kalman_row_state)
    {
      const auto &state = kalman_row_state;
      const double qop_t = state[TpcTrackKalmanFitter::QOverPt];
      const double omega = 0.003 * tracklet.kalman.bfield_t * qop_t;
      m_track.kalman_qop_t = static_cast<float>(qop_t);
      m_track.kalman_omega = static_cast<float>(omega);
      if (std::abs(omega) > 1.0e-12)
      {
        const double center_x = state[TpcTrackKalmanFitter::X] -
                                std::sin(state[TpcTrackKalmanFitter::Phi]) / omega;
        const double center_y = state[TpcTrackKalmanFitter::Y] +
                                std::cos(state[TpcTrackKalmanFitter::Phi]) / omega;
        m_track.kalman_cx = static_cast<float>(center_x);
        m_track.kalman_cy = static_cast<float>(center_y);
        m_track.kalman_radius = static_cast<float>(std::abs(1.0 / omega));
      }
    }
  }

  m_track.fit_chi2 = std::isfinite(tracklet.fit_chi2)
                         ? static_cast<float>(tracklet.fit_chi2)
                         : nan;
  m_track.fit_ndf = tracklet.fit_ndf;
  m_track.quality = std::isfinite(tracklet.fit_chi2_ndf)
                        ? static_cast<float>(tracklet.fit_chi2_ndf)
                        : nan;

  m_track.truth_px = static_cast<float>(tracklet.truth_momentum.x);
  m_track.truth_py = static_cast<float>(tracklet.truth_momentum.y);
  m_track.truth_pz = static_cast<float>(tracklet.truth_momentum.z);
  m_track.cos_mom_truth = static_cast<float>(vector_cosine(row_momentum, tracklet.truth_momentum));

  m_track_tree->Fill();
  ++m_counter_tracks_written;
}

void TpcV0CandidateTree::fill_cluster_residual_rows(const Tracklet &tracklet,
                                                    const Vec3 & /*primary_vertex*/,
                                                    const int run_number,
                                                    const int event_number)
{
  if (!m_cluster_residual_tree || tracklet.points.empty())
  {
    return;
  }

  const float nan = quiet_nan();
  double fit_chi2 = std::numeric_limits<double>::quiet_NaN();
  int fit_ndf = 0;

  const double sigma_rphi = std::max(m_kalman_config.meas_sigma_rphi_cm,
                                     m_kalman_config.min_measurement_sigma_cm);
  const double sigma_z = std::max(m_kalman_config.meas_sigma_z_cm,
                                  m_kalman_config.min_measurement_sigma_cm);

  auto helix_residual = [&](const TruthPoint &point,
                            double &previous_theta,
                            bool &have_previous_theta,
                            Vec3 &fit_position,
                            double &residual_r,
                            double &residual_rphi,
                            double &residual_z) -> bool
  {
    if (!tracklet.has_helix)
    {
      return false;
    }

    const Vec3 &cluster = point.position;
    const double cluster_r = pt(cluster);
    double theta_hint = std::atan2(cluster.y - tracklet.helix.cy,
                                   cluster.x - tracklet.helix.cx);
    theta_hint = unwrap_to_near(theta_hint,
                                have_previous_theta ? previous_theta
                                                    : tracklet.helix.theta_first);

    double theta = theta_hint;
    if (!helix_point_at_beam_radius(tracklet.helix, cluster_r, theta_hint, theta, fit_position))
    {
      fit_position = helix_point(tracklet.helix, theta_hint);
      theta = theta_hint;
      if (!finite(fit_position))
      {
        return false;
      }
    }

    previous_theta = theta;
    have_previous_theta = true;

    const double fit_r = pt(fit_position);
    const double cluster_phi = std::atan2(cluster.y, cluster.x);
    const double fit_phi = std::atan2(fit_position.y, fit_position.x);
    residual_r = cluster_r - fit_r;
    residual_rphi = cluster_r * normalize_phi(cluster_phi - fit_phi);
    residual_z = cluster.z - fit_position.z;
    return std::isfinite(residual_r) &&
           std::isfinite(residual_rphi) &&
           std::isfinite(residual_z);
  };

  if (tracklet.has_kalman)
  {
    fit_chi2 = tracklet.kalman.chi2;
    fit_ndf = tracklet.kalman.ndof;
  }
  else if (tracklet.has_helix)
  {
    double chi2 = 0.0;
    int nresiduals = 0;
    double previous_theta = tracklet.helix.theta_first;
    bool have_previous_theta = false;
    for (const auto &point : tracklet.points)
    {
      Vec3 fit_position;
      double residual_r = 0.0;
      double residual_rphi = 0.0;
      double residual_z = 0.0;
      if (!helix_residual(point, previous_theta, have_previous_theta,
                          fit_position, residual_r, residual_rphi, residual_z))
      {
        continue;
      }
      chi2 += square(residual_rphi / sigma_rphi) + square(residual_z / sigma_z);
      nresiduals += 2;
    }
    fit_chi2 = chi2;
    fit_ndf = std::max(0, nresiduals - 5);
  }

  double previous_theta = tracklet.has_helix ? tracklet.helix.theta_first : 0.0;
  bool have_previous_theta = false;
  for (std::size_t index = 0; index < tracklet.points.size(); ++index)
  {
    const auto &point = tracklet.points[index];
    const Vec3 &cluster = point.position;
    Vec3 fit_position{nan, nan, nan};
    double residual_r = std::numeric_limits<double>::quiet_NaN();
    double residual_rphi = std::numeric_limits<double>::quiet_NaN();
    double residual_z = std::numeric_limits<double>::quiet_NaN();

    if (tracklet.has_kalman && index < tracklet.kalman.states_smoothed.size())
    {
      fit_position = TpcTrackKalmanFitter::state_position(tracklet.kalman.states_smoothed[index]);
      const Vec3 delta = subtract(cluster, fit_position);
      const double fit_phi = std::atan2(fit_position.y, fit_position.x);
      residual_r = std::cos(fit_phi) * delta.x + std::sin(fit_phi) * delta.y;
      residual_rphi = -std::sin(fit_phi) * delta.x + std::cos(fit_phi) * delta.y;
      residual_z = delta.z;
    }
    else if (tracklet.has_helix)
    {
      helix_residual(point, previous_theta, have_previous_theta,
                     fit_position, residual_r, residual_rphi, residual_z);
    }

    reset_cluster_residual_row();
    m_cluster_residual.run = run_number;
    m_cluster_residual.evt = event_number;
    m_cluster_residual.track_id = tracklet.track_id;
    m_cluster_residual.charge = tracklet.charge;
    m_cluster_residual.side = (cluster.z >= 0.0) ? 1 : -1;
    m_cluster_residual.layer = point.layer;
    m_cluster_residual.cluster_index = static_cast<int>(index);
    m_cluster_residual.ntp_cluster = tracklet.npoints;
    m_cluster_residual.npoints = tracklet.npoints;
    m_cluster_residual.has_helix = tracklet.has_helix ? 1 : 0;
    m_cluster_residual.has_kalman = tracklet.has_kalman ? 1 : 0;

    m_cluster_residual.cluster_x = static_cast<float>(cluster.x);
    m_cluster_residual.cluster_y = static_cast<float>(cluster.y);
    m_cluster_residual.cluster_z = static_cast<float>(cluster.z);
    m_cluster_residual.cluster_r = static_cast<float>(pt(cluster));
    m_cluster_residual.cluster_phi = static_cast<float>(std::atan2(cluster.y, cluster.x));

    m_cluster_residual.fit_x = static_cast<float>(fit_position.x);
    m_cluster_residual.fit_y = static_cast<float>(fit_position.y);
    m_cluster_residual.fit_z = static_cast<float>(fit_position.z);
    m_cluster_residual.fit_r = static_cast<float>(pt(fit_position));
    m_cluster_residual.fit_phi = static_cast<float>(std::atan2(fit_position.y, fit_position.x));

    m_cluster_residual.residual_x = static_cast<float>(cluster.x - fit_position.x);
    m_cluster_residual.residual_y = static_cast<float>(cluster.y - fit_position.y);
    m_cluster_residual.residual_z = static_cast<float>(residual_z);
    m_cluster_residual.residual_r = static_cast<float>(residual_r);
    m_cluster_residual.residual_rphi = static_cast<float>(residual_rphi);

    m_cluster_residual.fit_chi2 = static_cast<float>(fit_chi2);
    m_cluster_residual.fit_ndf = fit_ndf;
    m_cluster_residual.fit_chi2_ndf =
        (fit_ndf > 0 && std::isfinite(fit_chi2)) ? static_cast<float>(fit_chi2 / fit_ndf) : nan;

    m_cluster_residual_tree->Fill();
    ++m_counter_cluster_residuals_written;
  }
}

void TpcV0CandidateTree::assign_fit_quality(Tracklet &tracklet) const
{
  tracklet.fit_chi2 = std::numeric_limits<double>::quiet_NaN();
  tracklet.fit_ndf = 0;
  tracklet.fit_chi2_ndf = std::numeric_limits<double>::quiet_NaN();

  if (tracklet.has_kalman)
  {
    tracklet.fit_chi2 = tracklet.kalman.chi2;
    tracklet.fit_ndf = tracklet.kalman.ndof;
  }
  else if (tracklet.has_helix)
  {
    const double sigma_rphi = std::max(m_kalman_config.meas_sigma_rphi_cm,
                                       m_kalman_config.min_measurement_sigma_cm);
    const double sigma_z = std::max(m_kalman_config.meas_sigma_z_cm,
                                    m_kalman_config.min_measurement_sigma_cm);

    double chi2 = 0.0;
    int nresiduals = 0;
    double previous_theta = tracklet.helix.theta_first;
    bool have_previous_theta = false;
    for (const auto &point : tracklet.points)
    {
      Vec3 fit_position;
      double residual_r = 0.0;
      double residual_rphi = 0.0;
      double residual_z = 0.0;
      if (!helix_cluster_residual(tracklet.helix, point, previous_theta, have_previous_theta,
                                  fit_position, residual_r, residual_rphi, residual_z))
      {
        continue;
      }
      chi2 += square(residual_rphi / sigma_rphi) + square(residual_z / sigma_z);
      nresiduals += 2;
    }

    tracklet.fit_chi2 = chi2;
    tracklet.fit_ndf = std::max(0, nresiduals - 5);
  }

  if (tracklet.fit_ndf > 0 && std::isfinite(tracklet.fit_chi2))
  {
    tracklet.fit_chi2_ndf = tracklet.fit_chi2 / tracklet.fit_ndf;
  }
}

void TpcV0CandidateTree::reset_pair_row()
{
  m_pair = {};
  const float nan = quiet_nan();
  m_pair.dca_xy1 = nan;
  m_pair.dca_z1 = nan;
  m_pair.dca_xy2 = nan;
  m_pair.dca_z2 = nan;
  m_pair.pairDCA = nan;
  m_pair.alpha = nan;
  m_pair.qT = nan;
  m_pair.dedx_1 = nan;
  m_pair.dedx_2 = nan;
  m_pair.cosThetaReco = nan;
  m_pair.Lproj = nan;
  m_pair.pca_x = nan;
  m_pair.pca_y = nan;
  m_pair.pca_z = nan;
  m_pair.pca1_x = nan;
  m_pair.pca1_y = nan;
  m_pair.pca1_z = nan;
  m_pair.pca2_x = nan;
  m_pair.pca2_y = nan;
  m_pair.pca2_z = nan;
  m_pair.v0_px = nan;
  m_pair.v0_py = nan;
  m_pair.v0_pz = nan;
  m_pair.v0_pt = nan;
  m_pair.mass_Kshort = nan;
  m_pair.mass_Lambda = nan;
  m_pair.mass_AntiLambda = nan;
  m_pair.true_decay_x = nan;
  m_pair.true_decay_y = nan;
  m_pair.true_decay_z = nan;
  m_pair.pca_to_true_3d = nan;
  m_pair.pca_to_true_xy = nan;
  m_pair.pca_to_true_z = nan;
  m_pair.truth_alpha = nan;
  m_pair.truth_qT = nan;
  m_pair.delta_alpha = nan;
  m_pair.delta_qT = nan;
  m_pair.truth_px1 = nan;
  m_pair.truth_py1 = nan;
  m_pair.truth_pz1 = nan;
  m_pair.truth_px2 = nan;
  m_pair.truth_py2 = nan;
  m_pair.truth_pz2 = nan;
  m_pair.cos_mom1_truth = nan;
  m_pair.cos_mom2_truth = nan;
  m_pair.pca_theta1 = nan;
  m_pair.pca_theta2 = nan;
  m_pair.kalman_chi2_1 = nan;
  m_pair.kalman_chi2_2 = nan;
  m_pair.kalman_chi2_ndf1 = nan;
  m_pair.kalman_chi2_ndf2 = nan;
  m_pair.quality1 = nan;
  m_pair.quality2 = nan;
  m_pair.kalman_ndof1 = 0;
  m_pair.kalman_ndof2 = 0;
}

void TpcV0CandidateTree::reset_track_row()
{
  m_track = {};
  const float nan = quiet_nan();
  m_track.px = nan;
  m_track.py = nan;
  m_track.pz = nan;
  m_track.pt = nan;
  m_track.p = nan;
  m_track.eta = nan;
  m_track.dedx = nan;
  m_track.x = nan;
  m_track.y = nan;
  m_track.z = nan;
  m_track.first_x = nan;
  m_track.first_y = nan;
  m_track.first_z = nan;
  m_track.first_r = nan;
  m_track.last_x = nan;
  m_track.last_y = nan;
  m_track.last_z = nan;
  m_track.last_r = nan;
  m_track.dca_xy = nan;
  m_track.dca_z = nan;
  m_track.vertex_x = nan;
  m_track.vertex_y = nan;
  m_track.vertex_z = nan;
  m_track.vertex_z_rms = nan;
  m_track.pca_x = nan;
  m_track.pca_y = nan;
  m_track.pca_z = nan;
  m_track.rDCA_zero = nan;
  m_track.zDCA = nan;
  m_track.helix_cx = nan;
  m_track.helix_cy = nan;
  m_track.helix_radius = nan;
  m_track.helix_z0 = nan;
  m_track.helix_pitch = nan;
  m_track.helix_theta_first = nan;
  m_track.helix_theta_last = nan;
  m_track.helix_direction = nan;
  m_track.helix_search_anchored = 0;
  m_track.helix_anchor_point_index = -1;
  m_track.helix_anchor_theta = nan;
  m_track.helix_anchor_path_cm = nan;
  m_track.helix_anchor_residual_cm = nan;
  m_track.helix_search_theta_min = nan;
  m_track.helix_search_theta_max = nan;
  m_track.helix_search_upstream_cm = nan;
  m_track.helix_search_downstream_cm = nan;
  m_track.kalman_chi2 = nan;
  m_track.kalman_ndof = 0;
  m_track.kalman_qop_t = nan;
  m_track.kalman_omega = nan;
  m_track.kalman_cx = nan;
  m_track.kalman_cy = nan;
  m_track.kalman_radius = nan;
  m_track.fit_chi2 = nan;
  m_track.fit_ndf = 0;
  m_track.quality = nan;
  m_track.truth_px = nan;
  m_track.truth_py = nan;
  m_track.truth_pz = nan;
  m_track.cos_mom_truth = nan;
}

void TpcV0CandidateTree::reset_cluster_residual_row()
{
  m_cluster_residual = {};
  const float nan = quiet_nan();
  m_cluster_residual.cluster_x = nan;
  m_cluster_residual.cluster_y = nan;
  m_cluster_residual.cluster_z = nan;
  m_cluster_residual.cluster_r = nan;
  m_cluster_residual.cluster_phi = nan;
  m_cluster_residual.fit_x = nan;
  m_cluster_residual.fit_y = nan;
  m_cluster_residual.fit_z = nan;
  m_cluster_residual.fit_r = nan;
  m_cluster_residual.fit_phi = nan;
  m_cluster_residual.residual_x = nan;
  m_cluster_residual.residual_y = nan;
  m_cluster_residual.residual_z = nan;
  m_cluster_residual.residual_r = nan;
  m_cluster_residual.residual_rphi = nan;
  m_cluster_residual.fit_chi2 = nan;
  m_cluster_residual.fit_ndf = 0;
  m_cluster_residual.fit_chi2_ndf = nan;
}

void TpcV0CandidateTree::create_branches()
{
  m_pair_tree->Branch("run", &m_pair.run, "run/I");
  m_pair_tree->Branch("evt", &m_pair.evt, "evt/I");
  m_pair_tree->Branch("cross1", &m_pair.cross1, "cross1/S");
  m_pair_tree->Branch("cross2", &m_pair.cross2, "cross2/S");
  m_pair_tree->Branch("px1", &m_pair.px1, "px1/F");
  m_pair_tree->Branch("py1", &m_pair.py1, "py1/F");
  m_pair_tree->Branch("pz1", &m_pair.pz1, "pz1/F");
  m_pair_tree->Branch("px2", &m_pair.px2, "px2/F");
  m_pair_tree->Branch("py2", &m_pair.py2, "py2/F");
  m_pair_tree->Branch("pz2", &m_pair.pz2, "pz2/F");
  m_pair_tree->Branch("dca_xy1", &m_pair.dca_xy1, "dca_xy1/F");
  m_pair_tree->Branch("dca_z1", &m_pair.dca_z1, "dca_z1/F");
  m_pair_tree->Branch("dca_xy2", &m_pair.dca_xy2, "dca_xy2/F");
  m_pair_tree->Branch("dca_z2", &m_pair.dca_z2, "dca_z2/F");
  m_pair_tree->Branch("pairDCA", &m_pair.pairDCA, "pairDCA/F");
  m_pair_tree->Branch("alpha", &m_pair.alpha, "alpha/F");
  m_pair_tree->Branch("qT", &m_pair.qT, "qT/F");
  m_pair_tree->Branch("charge1", &m_pair.charge1, "charge1/F");
  m_pair_tree->Branch("charge2", &m_pair.charge2, "charge2/F");
  m_pair_tree->Branch("dedx_1", &m_pair.dedx_1, "dedx_1/F");
  m_pair_tree->Branch("dedx_2", &m_pair.dedx_2, "dedx_2/F");
  m_pair_tree->Branch("cosThetaReco", &m_pair.cosThetaReco, "cosThetaReco/F");
  m_pair_tree->Branch("Lproj", &m_pair.Lproj, "Lproj/F");
  m_pair_tree->Branch("pca_x", &m_pair.pca_x, "pca_x/F");
  m_pair_tree->Branch("pca_y", &m_pair.pca_y, "pca_y/F");
  m_pair_tree->Branch("pca_z", &m_pair.pca_z, "pca_z/F");
  m_pair_tree->Branch("pca1_x", &m_pair.pca1_x, "pca1_x/F");
  m_pair_tree->Branch("pca1_y", &m_pair.pca1_y, "pca1_y/F");
  m_pair_tree->Branch("pca1_z", &m_pair.pca1_z, "pca1_z/F");
  m_pair_tree->Branch("pca2_x", &m_pair.pca2_x, "pca2_x/F");
  m_pair_tree->Branch("pca2_y", &m_pair.pca2_y, "pca2_y/F");
  m_pair_tree->Branch("pca2_z", &m_pair.pca2_z, "pca2_z/F");
  m_pair_tree->Branch("v0_px", &m_pair.v0_px, "v0_px/F");
  m_pair_tree->Branch("v0_py", &m_pair.v0_py, "v0_py/F");
  m_pair_tree->Branch("v0_pz", &m_pair.v0_pz, "v0_pz/F");
  m_pair_tree->Branch("v0_pt", &m_pair.v0_pt, "v0_pt/F");
  m_pair_tree->Branch("mass_Kshort", &m_pair.mass_Kshort, "mass_Kshort/F");
  m_pair_tree->Branch("mass_Lambda", &m_pair.mass_Lambda, "mass_Lambda/F");
  m_pair_tree->Branch("mass_AntiLambda", &m_pair.mass_AntiLambda, "mass_AntiLambda/F");
  m_pair_tree->Branch("true_decay_x", &m_pair.true_decay_x, "true_decay_x/F");
  m_pair_tree->Branch("true_decay_y", &m_pair.true_decay_y, "true_decay_y/F");
  m_pair_tree->Branch("true_decay_z", &m_pair.true_decay_z, "true_decay_z/F");
  m_pair_tree->Branch("pca_to_true_3d", &m_pair.pca_to_true_3d, "pca_to_true_3d/F");
  m_pair_tree->Branch("pca_to_true_xy", &m_pair.pca_to_true_xy, "pca_to_true_xy/F");
  m_pair_tree->Branch("pca_to_true_z", &m_pair.pca_to_true_z, "pca_to_true_z/F");
  m_pair_tree->Branch("truth_alpha", &m_pair.truth_alpha, "truth_alpha/F");
  m_pair_tree->Branch("truth_qT", &m_pair.truth_qT, "truth_qT/F");
  m_pair_tree->Branch("delta_alpha", &m_pair.delta_alpha, "delta_alpha/F");
  m_pair_tree->Branch("delta_qT", &m_pair.delta_qT, "delta_qT/F");
  m_pair_tree->Branch("truth_px1", &m_pair.truth_px1, "truth_px1/F");
  m_pair_tree->Branch("truth_py1", &m_pair.truth_py1, "truth_py1/F");
  m_pair_tree->Branch("truth_pz1", &m_pair.truth_pz1, "truth_pz1/F");
  m_pair_tree->Branch("truth_px2", &m_pair.truth_px2, "truth_px2/F");
  m_pair_tree->Branch("truth_py2", &m_pair.truth_py2, "truth_py2/F");
  m_pair_tree->Branch("truth_pz2", &m_pair.truth_pz2, "truth_pz2/F");
  m_pair_tree->Branch("cos_mom1_truth", &m_pair.cos_mom1_truth, "cos_mom1_truth/F");
  m_pair_tree->Branch("cos_mom2_truth", &m_pair.cos_mom2_truth, "cos_mom2_truth/F");
  m_pair_tree->Branch("pca_theta1", &m_pair.pca_theta1, "pca_theta1/F");
  m_pair_tree->Branch("pca_theta2", &m_pair.pca_theta2, "pca_theta2/F");
  m_pair_tree->Branch("kalman_chi2_1", &m_pair.kalman_chi2_1, "kalman_chi2_1/F");
  m_pair_tree->Branch("kalman_chi2_2", &m_pair.kalman_chi2_2, "kalman_chi2_2/F");
  m_pair_tree->Branch("kalman_chi2_ndf1", &m_pair.kalman_chi2_ndf1, "kalman_chi2_ndf1/F");
  m_pair_tree->Branch("kalman_chi2_ndf2", &m_pair.kalman_chi2_ndf2, "kalman_chi2_ndf2/F");
  m_pair_tree->Branch("quality1", &m_pair.quality1, "quality1/F");
  m_pair_tree->Branch("quality2", &m_pair.quality2, "quality2/F");
  m_pair_tree->Branch("track_id1", &m_pair.track_id1, "track_id1/I");
  m_pair_tree->Branch("track_id2", &m_pair.track_id2, "track_id2/I");
  m_pair_tree->Branch("pid1", &m_pair.pid1, "pid1/I");
  m_pair_tree->Branch("pid2", &m_pair.pid2, "pid2/I");
  m_pair_tree->Branch("parent_id1", &m_pair.parent_id1, "parent_id1/I");
  m_pair_tree->Branch("parent_id2", &m_pair.parent_id2, "parent_id2/I");
  m_pair_tree->Branch("parent_pid", &m_pair.parent_pid, "parent_pid/I");
  m_pair_tree->Branch("kalman_ndof1", &m_pair.kalman_ndof1, "kalman_ndof1/I");
  m_pair_tree->Branch("kalman_ndof2", &m_pair.kalman_ndof2, "kalman_ndof2/I");
  m_pair_tree->Branch("npoints1", &m_pair.npoints1, "npoints1/S");
  m_pair_tree->Branch("npoints2", &m_pair.npoints2, "npoints2/S");

  m_track_tree->Branch("run", &m_track.run, "run/I");
  m_track_tree->Branch("evt", &m_track.evt, "evt/I");
  m_track_tree->Branch("track_id", &m_track.track_id, "track_id/I");
  m_track_tree->Branch("shower_id", &m_track.shower_id, "shower_id/I");
  m_track_tree->Branch("pid", &m_track.pid, "pid/I");
  m_track_tree->Branch("parent_id", &m_track.parent_id, "parent_id/I");
  m_track_tree->Branch("parent_pid", &m_track.parent_pid, "parent_pid/I");
  m_track_tree->Branch("charge", &m_track.charge, "charge/D");
  m_track_tree->Branch("side", &m_track.side, "side/I");
  m_track_tree->Branch("npoints", &m_track.npoints, "npoints/I");
  m_track_tree->Branch("ntpc_clusters", &m_track.ntpc_clusters, "ntpc_clusters/i");
  m_track_tree->Branch("has_helix", &m_track.has_helix, "has_helix/I");
  m_track_tree->Branch("has_kalman", &m_track.has_kalman, "has_kalman/I");
  m_track_tree->Branch("is_primary", &m_track.is_primary, "is_primary/I");
  m_track_tree->Branch("px", &m_track.px, "px/D");
  m_track_tree->Branch("py", &m_track.py, "py/D");
  m_track_tree->Branch("pz", &m_track.pz, "pz/D");
  m_track_tree->Branch("pt", &m_track.pt, "pt/D");
  m_track_tree->Branch("p", &m_track.p, "p/D");
  m_track_tree->Branch("eta", &m_track.eta, "eta/D");
  m_track_tree->Branch("dedx", &m_track.dedx, "dedx/D");
  m_track_tree->Branch("x", &m_track.x, "x/F");
  m_track_tree->Branch("y", &m_track.y, "y/F");
  m_track_tree->Branch("z", &m_track.z, "z/F");
  m_track_tree->Branch("first_x", &m_track.first_x, "first_x/F");
  m_track_tree->Branch("first_y", &m_track.first_y, "first_y/F");
  m_track_tree->Branch("first_z", &m_track.first_z, "first_z/F");
  m_track_tree->Branch("first_r", &m_track.first_r, "first_r/F");
  m_track_tree->Branch("last_x", &m_track.last_x, "last_x/F");
  m_track_tree->Branch("last_y", &m_track.last_y, "last_y/F");
  m_track_tree->Branch("last_z", &m_track.last_z, "last_z/F");
  m_track_tree->Branch("last_r", &m_track.last_r, "last_r/F");
  m_track_tree->Branch("dca_xy", &m_track.dca_xy, "dca_xy/F");
  m_track_tree->Branch("dca_z", &m_track.dca_z, "dca_z/F");
  m_track_tree->Branch("vertex_x", &m_track.vertex_x, "vertex_x/D");
  m_track_tree->Branch("vertex_y", &m_track.vertex_y, "vertex_y/D");
  m_track_tree->Branch("vertex_z", &m_track.vertex_z, "vertex_z/D");
  m_track_tree->Branch("vertex_from_upstream", &m_track.vertex_from_upstream,
                       "vertex_from_upstream/I");
  m_track_tree->Branch("vertex_z_rms", &m_track.vertex_z_rms, "vertex_z_rms/D");
  m_track_tree->Branch("vertex_ntracks", &m_track.vertex_ntracks, "vertex_ntracks/i");
  m_track_tree->Branch("pca_x", &m_track.pca_x, "pca_x/D");
  m_track_tree->Branch("pca_y", &m_track.pca_y, "pca_y/D");
  m_track_tree->Branch("pca_z", &m_track.pca_z, "pca_z/D");
  m_track_tree->Branch("rDCA_zero", &m_track.rDCA_zero, "rDCA_zero/D");
  m_track_tree->Branch("zDCA", &m_track.zDCA, "zDCA/D");
  m_track_tree->Branch("helix_cx", &m_track.helix_cx, "helix_cx/F");
  m_track_tree->Branch("helix_cy", &m_track.helix_cy, "helix_cy/F");
  m_track_tree->Branch("helix_radius", &m_track.helix_radius, "helix_radius/F");
  m_track_tree->Branch("helix_z0", &m_track.helix_z0, "helix_z0/F");
  m_track_tree->Branch("helix_pitch", &m_track.helix_pitch, "helix_pitch/F");
  m_track_tree->Branch("helix_theta_first", &m_track.helix_theta_first, "helix_theta_first/F");
  m_track_tree->Branch("helix_theta_last", &m_track.helix_theta_last, "helix_theta_last/F");
  m_track_tree->Branch("helix_direction", &m_track.helix_direction, "helix_direction/F");
  m_track_tree->Branch("helix_search_anchored", &m_track.helix_search_anchored,
                       "helix_search_anchored/I");
  m_track_tree->Branch("helix_anchor_point_index", &m_track.helix_anchor_point_index,
                       "helix_anchor_point_index/I");
  m_track_tree->Branch("helix_anchor_theta", &m_track.helix_anchor_theta,
                       "helix_anchor_theta/F");
  m_track_tree->Branch("helix_anchor_path_cm", &m_track.helix_anchor_path_cm,
                       "helix_anchor_path_cm/F");
  m_track_tree->Branch("helix_anchor_residual_cm", &m_track.helix_anchor_residual_cm,
                       "helix_anchor_residual_cm/F");
  m_track_tree->Branch("helix_search_theta_min", &m_track.helix_search_theta_min,
                       "helix_search_theta_min/F");
  m_track_tree->Branch("helix_search_theta_max", &m_track.helix_search_theta_max,
                       "helix_search_theta_max/F");
  m_track_tree->Branch("helix_search_upstream_cm", &m_track.helix_search_upstream_cm,
                       "helix_search_upstream_cm/F");
  m_track_tree->Branch("helix_search_downstream_cm", &m_track.helix_search_downstream_cm,
                       "helix_search_downstream_cm/F");
  m_track_tree->Branch("kalman_chi2", &m_track.kalman_chi2, "kalman_chi2/F");
  m_track_tree->Branch("kalman_ndof", &m_track.kalman_ndof, "kalman_ndof/I");
  m_track_tree->Branch("kalman_naccepted", &m_track.kalman_naccepted, "kalman_naccepted/i");
  m_track_tree->Branch("kalman_nrejected", &m_track.kalman_nrejected, "kalman_nrejected/i");
  m_track_tree->Branch("kalman_measurement_chi2", &m_track.kalman_measurement_chi2);
  m_track_tree->Branch("kalman_measurement_used", &m_track.kalman_measurement_used);
  if (m_kalman_config.collect_innovation_components)
  {
    m_track_tree->Branch("kalman_measurement_sigma_r", &m_track.kalman_measurement_sigma_r,
                         "kalman_measurement_sigma_r/F");
    m_track_tree->Branch("kalman_measurement_sigma_rphi", &m_track.kalman_measurement_sigma_rphi,
                         "kalman_measurement_sigma_rphi/F");
    m_track_tree->Branch("kalman_measurement_sigma_z", &m_track.kalman_measurement_sigma_z,
                         "kalman_measurement_sigma_z/F");
    m_track_tree->Branch("kalman_measurement_in_seed", &m_track.kalman_measurement_in_seed);
    m_track_tree->Branch("kalman_innovation_residual_r", &m_track.kalman_innovation_residual_r);
    m_track_tree->Branch("kalman_innovation_residual_rphi", &m_track.kalman_innovation_residual_rphi);
    m_track_tree->Branch("kalman_innovation_residual_z", &m_track.kalman_innovation_residual_z);
    m_track_tree->Branch("kalman_prediction_sigma_r", &m_track.kalman_prediction_sigma_r);
    m_track_tree->Branch("kalman_prediction_sigma_rphi", &m_track.kalman_prediction_sigma_rphi);
    m_track_tree->Branch("kalman_prediction_sigma_z", &m_track.kalman_prediction_sigma_z);
    m_track_tree->Branch("kalman_innovation_sigma_r", &m_track.kalman_innovation_sigma_r);
    m_track_tree->Branch("kalman_innovation_sigma_rphi", &m_track.kalman_innovation_sigma_rphi);
    m_track_tree->Branch("kalman_innovation_sigma_z", &m_track.kalman_innovation_sigma_z);
    m_track_tree->Branch("kalman_innovation_rho_r_rphi", &m_track.kalman_innovation_rho_r_rphi);
    m_track_tree->Branch("kalman_innovation_rho_r_z", &m_track.kalman_innovation_rho_r_z);
    m_track_tree->Branch("kalman_innovation_rho_rphi_z", &m_track.kalman_innovation_rho_rphi_z);
    m_track_tree->Branch("kalman_innovation_whitened_0", &m_track.kalman_innovation_whitened_0);
    m_track_tree->Branch("kalman_innovation_whitened_1", &m_track.kalman_innovation_whitened_1);
    m_track_tree->Branch("kalman_innovation_whitened_2", &m_track.kalman_innovation_whitened_2);
  }
  m_track_tree->Branch("kalman_qop_t", &m_track.kalman_qop_t, "kalman_qop_t/F");
  m_track_tree->Branch("kalman_omega", &m_track.kalman_omega, "kalman_omega/F");
  m_track_tree->Branch("kalman_cx", &m_track.kalman_cx, "kalman_cx/F");
  m_track_tree->Branch("kalman_cy", &m_track.kalman_cy, "kalman_cy/F");
  m_track_tree->Branch("kalman_radius", &m_track.kalman_radius, "kalman_radius/F");
  m_track_tree->Branch("fit_chi2", &m_track.fit_chi2, "fit_chi2/F");
  m_track_tree->Branch("fit_ndf", &m_track.fit_ndf, "fit_ndf/I");
  m_track_tree->Branch("quality", &m_track.quality, "quality/F");
  m_track_tree->Branch("truth_px", &m_track.truth_px, "truth_px/F");
  m_track_tree->Branch("truth_py", &m_track.truth_py, "truth_py/F");
  m_track_tree->Branch("truth_pz", &m_track.truth_pz, "truth_pz/F");
  m_track_tree->Branch("cos_mom_truth", &m_track.cos_mom_truth, "cos_mom_truth/F");
  m_track_tree->Branch("cluster_index", &m_track.cluster_index);
  m_track_tree->Branch("cluster_side", &m_track.cluster_side);
  m_track_tree->Branch("layer", &m_track.layer);
  m_track_tree->Branch("cluster_z", &m_track.cluster_z);
  m_track_tree->Branch("cluster_r", &m_track.cluster_r);
  m_track_tree->Branch("cluster_phi", &m_track.cluster_phi);
  m_track_tree->Branch("residual_z", &m_track.residual_z);
  m_track_tree->Branch("residual_r", &m_track.residual_r);
  m_track_tree->Branch("residual_rphi", &m_track.residual_rphi);

  if (!m_cluster_residual_tree)
  {
    return;
  }

  m_cluster_residual_tree->Branch("run", &m_cluster_residual.run, "run/I");
  m_cluster_residual_tree->Branch("evt", &m_cluster_residual.evt, "evt/I");
  m_cluster_residual_tree->Branch("track_id", &m_cluster_residual.track_id, "track_id/I");
  m_cluster_residual_tree->Branch("charge", &m_cluster_residual.charge, "charge/I");
  m_cluster_residual_tree->Branch("side", &m_cluster_residual.side, "side/I");
  m_cluster_residual_tree->Branch("layer", &m_cluster_residual.layer, "layer/I");
  m_cluster_residual_tree->Branch("cluster_index", &m_cluster_residual.cluster_index, "cluster_index/I");
  m_cluster_residual_tree->Branch("ntp_cluster", &m_cluster_residual.ntp_cluster, "ntp_cluster/I");
  m_cluster_residual_tree->Branch("npoints", &m_cluster_residual.npoints, "npoints/I");
  m_cluster_residual_tree->Branch("has_helix", &m_cluster_residual.has_helix, "has_helix/I");
  m_cluster_residual_tree->Branch("has_kalman", &m_cluster_residual.has_kalman, "has_kalman/I");
  m_cluster_residual_tree->Branch("cluster_x", &m_cluster_residual.cluster_x, "cluster_x/F");
  m_cluster_residual_tree->Branch("cluster_y", &m_cluster_residual.cluster_y, "cluster_y/F");
  m_cluster_residual_tree->Branch("cluster_z", &m_cluster_residual.cluster_z, "cluster_z/F");
  m_cluster_residual_tree->Branch("cluster_r", &m_cluster_residual.cluster_r, "cluster_r/F");
  m_cluster_residual_tree->Branch("cluster_phi", &m_cluster_residual.cluster_phi, "cluster_phi/F");
  m_cluster_residual_tree->Branch("fit_x", &m_cluster_residual.fit_x, "fit_x/F");
  m_cluster_residual_tree->Branch("fit_y", &m_cluster_residual.fit_y, "fit_y/F");
  m_cluster_residual_tree->Branch("fit_z", &m_cluster_residual.fit_z, "fit_z/F");
  m_cluster_residual_tree->Branch("fit_r", &m_cluster_residual.fit_r, "fit_r/F");
  m_cluster_residual_tree->Branch("fit_phi", &m_cluster_residual.fit_phi, "fit_phi/F");
  m_cluster_residual_tree->Branch("residual_x", &m_cluster_residual.residual_x, "residual_x/F");
  m_cluster_residual_tree->Branch("residual_y", &m_cluster_residual.residual_y, "residual_y/F");
  m_cluster_residual_tree->Branch("residual_z", &m_cluster_residual.residual_z, "residual_z/F");
  m_cluster_residual_tree->Branch("residual_r", &m_cluster_residual.residual_r, "residual_r/F");
  m_cluster_residual_tree->Branch("residual_rphi", &m_cluster_residual.residual_rphi, "residual_rphi/F");
  m_cluster_residual_tree->Branch("fit_chi2", &m_cluster_residual.fit_chi2, "fit_chi2/F");
  m_cluster_residual_tree->Branch("fit_ndf", &m_cluster_residual.fit_ndf, "fit_ndf/I");
  m_cluster_residual_tree->Branch("fit_chi2_ndf", &m_cluster_residual.fit_chi2_ndf, "fit_chi2_ndf/F");
}

int TpcV0CandidateTree::pdg_charge(const int pid)
{
  const int apid = std::abs(pid);
  int charge = 0;
  switch (apid)
  {
  case 11:
  case 13:
    charge = -1;
    break;
  case 211:
  case 321:
  case 2212:
  case 3222:
    charge = 1;
    break;
  case 3112:
  case 3312:
  case 3334:
    charge = -1;
    break;
  default:
    charge = 0;
    break;
  }
  return (pid < 0) ? -charge : charge;
}

float TpcV0CandidateTree::quiet_nan()
{
  return std::numeric_limits<float>::quiet_NaN();
}

bool TpcV0CandidateTree::parse_point_order(const std::string &mode, PointOrder &order)
{
  return TpcTrackHelixFitter::parse_point_order(mode, order);
}

bool TpcV0CandidateTree::finite(const Vec3 &value)
{
  return TpcTrackHelixFitter::finite(value);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::add(const Vec3 &lhs, const Vec3 &rhs)
{
  return TpcTrackHelixFitter::add(lhs, rhs);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::subtract(const Vec3 &lhs, const Vec3 &rhs)
{
  return TpcTrackHelixFitter::subtract(lhs, rhs);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::scale(const Vec3 &value, const double factor)
{
  return TpcTrackHelixFitter::scale(value, factor);
}

double TpcV0CandidateTree::dot(const Vec3 &lhs, const Vec3 &rhs)
{
  return TpcTrackHelixFitter::dot(lhs, rhs);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::cross(const Vec3 &lhs, const Vec3 &rhs)
{
  return {
      lhs.y * rhs.z - lhs.z * rhs.y,
      lhs.z * rhs.x - lhs.x * rhs.z,
      lhs.x * rhs.y - lhs.y * rhs.x};
}

double TpcV0CandidateTree::norm(const Vec3 &value)
{
  return TpcTrackHelixFitter::norm(value);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::unit(const Vec3 &value)
{
  return TpcTrackHelixFitter::unit(value);
}

double TpcV0CandidateTree::pt(const Vec3 &value)
{
  return TpcTrackHelixFitter::pt(value);
}

double TpcV0CandidateTree::distance(const Vec3 &lhs, const Vec3 &rhs)
{
  return TpcTrackHelixFitter::distance(lhs, rhs);
}

double TpcV0CandidateTree::vector_cosine(const Vec3 &lhs, const Vec3 &rhs)
{
  return TpcTrackHelixFitter::vector_cosine(lhs, rhs);
}

bool TpcV0CandidateTree::fit_circle_least_squares(const std::vector<TruthPoint> &points,
                                                  const std::size_t nfit,
                                                  double &cx,
                                                  double &cy,
                                                  double &radius)
{
  return TpcTrackHelixFitter::fit_circle_least_squares(points, nfit, cx, cy, radius);
}

void TpcV0CandidateTree::order_track_points(std::vector<TruthPoint> &points, const PointOrder order)
{
  TpcTrackHelixFitter::order_points(points, order);
}

bool TpcV0CandidateTree::fit_helix(const std::vector<TruthPoint> &points,
                                   const int fit_first_points,
                                   const int charge,
                                   const double bfield_t,
                                   HelixFit &helix)
{
  return TpcTrackHelixFitter::fit(points, fit_first_points, bfield_t, helix) &&
         TpcTrackHelixFitter::orient_to_charge(helix, charge);
}

bool TpcV0CandidateTree::fit_kalman(const std::vector<TruthPoint> &points,
                                    const int charge,
                                    TpcKalmanResult &kalman) const
{
  const auto start = std::chrono::steady_clock::now();
  const bool success = TpcTrackKalmanFitter::fit(
      points, charge, m_kalman_config, kalman, kPionMass);
  const double fit_seconds = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - start)
                                 .count();
  m_timing_kalman_fit_seconds += fit_seconds;
  m_timing_rkn_seconds += kalman.rkn_seconds;
  m_timing_rkn_propagations += kalman.rkn_propagations;
  m_timing_rkn_accepted_steps += kalman.rkn_accepted_steps;
  m_timing_rkn_rejected_trials += kalman.rkn_rejected_trials;
  m_timing_rkn_failures += kalman.rkn_failures;
  ++m_timing_kalman_fits;
  if (m_print_timing &&
      (m_timing_kalman_fits <= 3 || m_timing_kalman_fits % 10 == 0 || fit_seconds > 1.0))
  {
    std::cout << "[V0TimingFit] fit=" << m_timing_kalman_fits
              << " points=" << points.size()
              << " success=" << success
              << " fit_s=" << fit_seconds
              << " rkn_s=" << kalman.rkn_seconds
              << " rkn_propagations=" << kalman.rkn_propagations
              << " rkn_steps=" << kalman.rkn_accepted_steps
              << " rkn_retries=" << kalman.rkn_rejected_trials
              << " rkn_failures=" << kalman.rkn_failures
              << std::endl;
  }
  return success;
}

bool TpcV0CandidateTree::helix_from_state(const Vec3 &position, const Vec3 &momentum,
                                          const int charge, const double bfield_t,
                                          HelixFit &helix)
{
  return TpcTrackHelixFitter::from_state(position, momentum, charge, bfield_t, helix);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::helix_point(const HelixFit &helix, const double theta)
{
  return TpcTrackHelixFitter::point(helix, theta);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::helix_tangent(const HelixFit &helix, const double theta)
{
  return TpcTrackHelixFitter::tangent(helix, theta);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::helix_momentum(const HelixFit &helix, const double theta)
{
  return TpcTrackHelixFitter::momentum(helix, theta);
}

std::pair<double, double> TpcV0CandidateTree::theta_search_range(const HelixFit &helix,
                                                                 const double theta_extension,
                                                                 const double downstream_margin)
{
  return TpcTrackHelixFitter::theta_search_range(helix, theta_extension, downstream_margin);
}

bool TpcV0CandidateTree::line_line_pca(const Vec3 &pos1, const Vec3 &dir1,
                                       const Vec3 &pos2, const Vec3 &dir2,
                                       LinePca &pca, const bool normalize_dirs)
{
  return TpcTrackHelixFitter::line_line_pca(pos1, dir1, pos2, dir2, pca, normalize_dirs);
}

TpcV0CandidateTree::HelixPca TpcV0CandidateTree::refine_helix_pair(
    const HelixFit &helix1, const HelixFit &helix2,
    double theta1, double theta2,
    const double min1, const double max1,
    const double min2, const double max2,
    double max_step)
{
  return TpcTrackHelixFitter::refine_pair(helix1, helix2, theta1, theta2,
                                          min1, max1, min2, max2, max_step);
}

std::vector<TpcV0CandidateTree::HelixPca> TpcV0CandidateTree::helix_helix_pca_candidates(
    const HelixFit &helix1, const HelixFit &helix2,
    const double theta_extension,
    const int coarse_steps,
    const double downstream_margin,
    const int max_candidates)
{
  return TpcTrackHelixFitter::pca_candidates(helix1, helix2, theta_extension,
                                             coarse_steps, downstream_margin,
                                             max_candidates);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::kalman_point(const TpcKalmanResult &kalman,
                                                          const double s_cm,
                                                          const TpcKalmanConfig &config,
                                                          const Vec3 &reference_vertex)
{
  if (!kalman.success || kalman.states_smoothed.empty())
  {
    const double nan = quiet_nan();
    return {nan, nan, nan};
  }

  const auto state = TpcTrackKalmanFitter::propagate_state(
      TpcTrackKalmanFitter::propagation_state(kalman, reference_vertex), s_cm, config, kalman.mass_gev);
  return TpcTrackKalmanFitter::state_position(state);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::kalman_tangent(const TpcKalmanResult &kalman,
                                                            const double s_cm,
                                                            const TpcKalmanConfig &config,
                                                            const Vec3 &reference_vertex)
{
  if (!kalman.success || kalman.states_smoothed.empty())
  {
    const double nan = quiet_nan();
    return {nan, nan, nan};
  }

  const auto state = TpcTrackKalmanFitter::propagate_state(
      TpcTrackKalmanFitter::propagation_state(kalman, reference_vertex), s_cm, config, kalman.mass_gev);
  return TpcTrackKalmanFitter::state_tangent(state);
}

TpcV0CandidateTree::Vec3 TpcV0CandidateTree::kalman_momentum(const TpcKalmanResult &kalman,
                                                             const double s_cm,
                                                             const TpcKalmanConfig &config,
                                                             const Vec3 &reference_vertex)
{
  if (!kalman.success || kalman.states_smoothed.empty())
  {
    const double nan = quiet_nan();
    return {nan, nan, nan};
  }

  const auto state = TpcTrackKalmanFitter::propagate_state(
      TpcTrackKalmanFitter::propagation_state(kalman, reference_vertex), s_cm, config, kalman.mass_gev);
  return TpcTrackKalmanFitter::state_momentum(state);
}

TpcV0CandidateTree::KalmanPca TpcV0CandidateTree::refine_kalman_pair(
    const TpcKalmanResult &kalman1,
    const TpcKalmanResult &kalman2,
    const TpcKalmanConfig &config,
    const Vec3 &reference_vertex,
    double s1, double s2,
    const double min1, const double max1,
    const double min2, const double max2,
    double max_step,
    const int max_iterations)
{
  KalmanPca best;
  best.s1 = s1;
  best.s2 = s2;
  best.pca1 = kalman_point(kalman1, s1, config, reference_vertex);
  best.pca2 = kalman_point(kalman2, s2, config, reference_vertex);
  double best_dca2 = square(distance(best.pca1, best.pca2));

  for (int iter = 0; iter < std::max(1, max_iterations); ++iter)
  {
    const Vec3 pos1 = kalman_point(kalman1, s1, config, reference_vertex);
    const Vec3 pos2 = kalman_point(kalman2, s2, config, reference_vertex);
    const Vec3 tan1 = kalman_tangent(kalman1, s1, config, reference_vertex);
    const Vec3 tan2 = kalman_tangent(kalman2, s2, config, reference_vertex);

    LinePca line_pca;
    if (!line_line_pca(pos1, tan1, pos2, tan2, line_pca, false))
    {
      break;
    }

    const double step1 = std::clamp(line_pca.step1, -max_step, max_step);
    const double step2 = std::clamp(line_pca.step2, -max_step, max_step);
    if (std::abs(step1) < 1.0e-4 && std::abs(step2) < 1.0e-4)
    {
      break;
    }

    const double candidate_s1 = std::clamp(s1 + step1, min1, max1);
    const double candidate_s2 = std::clamp(s2 + step2, min2, max2);
    const Vec3 candidate_pos1 = kalman_point(kalman1, candidate_s1, config, reference_vertex);
    const Vec3 candidate_pos2 = kalman_point(kalman2, candidate_s2, config, reference_vertex);
    const double candidate_dca2 = square(distance(candidate_pos1, candidate_pos2));
    if (candidate_dca2 < best_dca2)
    {
      s1 = candidate_s1;
      s2 = candidate_s2;
      best_dca2 = candidate_dca2;
    }
    else
    {
      max_step *= 0.5;
      if (max_step < 1.0e-3)
      {
        break;
      }
    }
  }

  best.s1 = s1;
  best.s2 = s2;
  best.pca1 = kalman_point(kalman1, s1, config, reference_vertex);
  best.pca2 = kalman_point(kalman2, s2, config, reference_vertex);
  best.dca = distance(best.pca1, best.pca2);
  return best;
}

std::vector<TpcV0CandidateTree::KalmanPca> TpcV0CandidateTree::kalman_pca_candidates(
    const TpcKalmanResult &kalman1,
    const TpcKalmanResult &kalman2,
    const TpcKalmanConfig &config,
    const Vec3 &reference_vertex,
    const double max_upstream_cm,
    const double downstream_margin_cm,
    const int coarse_steps,
    const int max_candidates)
{
  std::vector<KalmanPca> candidates;
  if (!kalman1.success || !kalman2.success ||
      kalman1.states_smoothed.empty() || kalman2.states_smoothed.empty())
  {
    return candidates;
  }

  const double min1 = -std::abs(max_upstream_cm);
  const double max1 = std::abs(downstream_margin_cm);
  const double min2 = -std::abs(max_upstream_cm);
  const double max2 = std::abs(downstream_margin_cm);
  const int nsteps = std::max(8, coarse_steps);
  const int ncandidates = std::max(1, max_candidates);
  const bool use_fast_field_pca =
      config.magnetic_field != nullptr && config.rkn_fast_field_pca;
  TpcKalmanConfig search_config = config;
  if (use_fast_field_pca)
  {
    // Use a cheap uniform-Bz trajectory only to locate promising path lengths.
    // The selected candidates are refined below with the full field map.
    search_config.magnetic_field = nullptr;
    search_config.rkn_step_tolerance = 0.0;
    search_config.rkn_max_step_cm = std::max(20.0, std::abs(config.rkn_max_step_cm));
  }
  const TpcKalmanConfig &coarse_config = use_fast_field_pca ? search_config : config;

  std::vector<std::tuple<double, double, double>> seeds;
  const auto nsteps_size = static_cast<std::size_t>(nsteps);
  seeds.reserve(nsteps_size * nsteps_size);
  std::vector<double> s1_values(static_cast<std::size_t>(nsteps));
  std::vector<double> s2_values(static_cast<std::size_t>(nsteps));
  std::vector<Vec3> points1(static_cast<std::size_t>(nsteps));
  std::vector<Vec3> points2(static_cast<std::size_t>(nsteps));
  for (int i = 0; i < nsteps; ++i)
  {
    s1_values[static_cast<std::size_t>(i)] =
        min1 + (max1 - min1) * static_cast<double>(i) / static_cast<double>(nsteps - 1);
    s2_values[static_cast<std::size_t>(i)] =
        min2 + (max2 - min2) * static_cast<double>(i) / static_cast<double>(nsteps - 1);
    points1[static_cast<std::size_t>(i)] =
        kalman_point(kalman1, s1_values[static_cast<std::size_t>(i)], coarse_config, reference_vertex);
    points2[static_cast<std::size_t>(i)] =
        kalman_point(kalman2, s2_values[static_cast<std::size_t>(i)], coarse_config, reference_vertex);
  }

  for (int i = 0; i < nsteps; ++i)
  {
    const double s1 = s1_values[static_cast<std::size_t>(i)];
    const Vec3 &pos1 = points1[static_cast<std::size_t>(i)];
    for (int j = 0; j < nsteps; ++j)
    {
      const double s2 = s2_values[static_cast<std::size_t>(j)];
      const Vec3 &pos2 = points2[static_cast<std::size_t>(j)];
      const double d2 = square(distance(pos1, pos2));
      if (std::isfinite(d2))
      {
        seeds.emplace_back(d2, s1, s2);
      }
    }
  }

  std::sort(seeds.begin(), seeds.end(),
            [](const auto &lhs, const auto &rhs)
            { return std::get<0>(lhs) < std::get<0>(rhs); });

  const double max_step = std::max(max1 - min1, max2 - min2) / static_cast<double>(nsteps);
  const int nrefine = std::min(ncandidates, static_cast<int>(seeds.size()));
  candidates.reserve(static_cast<std::size_t>(nrefine));
  for (int index = 0; index < nrefine; ++index)
  {
    double seed_s1 = std::get<1>(seeds[static_cast<std::size_t>(index)]);
    double seed_s2 = std::get<2>(seeds[static_cast<std::size_t>(index)]);
    if (use_fast_field_pca)
    {
      const auto surrogate = refine_kalman_pair(
          kalman1, kalman2, search_config, reference_vertex,
          seed_s1, seed_s2, min1, max1, min2, max2, max_step, 12);
      seed_s1 = surrogate.s1;
      seed_s2 = surrogate.s2;
    }

    candidates.push_back(refine_kalman_pair(
        kalman1, kalman2, config, reference_vertex,
        seed_s1, seed_s2, min1, max1, min2, max2, max_step,
        use_fast_field_pca
            ? std::max(1, config.rkn_field_pca_refine_iterations)
            : 30));
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const KalmanPca &lhs, const KalmanPca &rhs)
            { return lhs.dca < rhs.dca; });
  return candidates;
}

std::pair<double, double> TpcV0CandidateTree::track_dca_to_vertex(const Vec3 &pos,
                                                                  const Vec3 &mom,
                                                                  const Vec3 &vertex)
{
  return TpcTrackHelixFitter::line_dca_to_vertex(pos, mom, vertex);
}

std::pair<double, double> TpcV0CandidateTree::helix_dca_to_vertex(const HelixFit &helix,
                                                                  const Vec3 &vertex)
{
  return TpcTrackHelixFitter::helix_dca_to_vertex(helix, vertex);
}

std::pair<double, double> TpcV0CandidateTree::fitted_track_dca_to_vertex(
    const Tracklet &tracklet,
    const Vec3 &vertex) const
{
  if (tracklet.has_kalman)
  {
    return TpcTrackKalmanFitter::dca_to_vertex(
        tracklet.kalman, vertex, &m_kalman_config);
  }
  if (tracklet.has_helix)
  {
    return helix_dca_to_vertex(tracklet.helix, vertex);
  }
  return track_dca_to_vertex(tracklet.position, tracklet.momentum, vertex);
}

bool TpcV0CandidateTree::armenteros(const Vec3 &pplus, const Vec3 &pminus,
                                    double &alpha, double &qt)
{
  const Vec3 v0p = add(pplus, pminus);
  const Vec3 direction = unit(v0p);
  if (!finite(direction))
  {
    return false;
  }

  const double pl_plus = dot(pplus, direction);
  const double pl_minus = dot(pminus, direction);
  const double denom = pl_plus + pl_minus;
  if (std::abs(denom) < 1e-10)
  {
    return false;
  }

  alpha = (pl_plus - pl_minus) / denom;
  qt = norm(subtract(pplus, scale(direction, pl_plus)));
  return true;
}

double TpcV0CandidateTree::invariant_mass(const Vec3 &mom1, const double mass1,
                                          const Vec3 &mom2, const double mass2)
{
  const double e1 = std::sqrt(dot(mom1, mom1) + square(mass1));
  const double e2 = std::sqrt(dot(mom2, mom2) + square(mass2));
  const Vec3 total_mom = add(mom1, mom2);
  const double mass2_total = square(e1 + e2) - dot(total_mom, total_mom);
  return (mass2_total > 0.0) ? std::sqrt(mass2_total) : 0.0;
}

bool TpcV0CandidateTree::passes_preselection(const Tracklet &track1, const Tracklet &track2,
                                             const Vec3 &primary_vertex) const
{
  if (m_pre_track_npoints_min > 0 &&
      (track1.npoints < m_pre_track_npoints_min || track2.npoints < m_pre_track_npoints_min))
  {
    return false;
  }

  if (m_pre_track_quality_max >= 0.0 &&
      (!std::isfinite(track1.fit_chi2_ndf) || !std::isfinite(track2.fit_chi2_ndf) ||
       track1.fit_chi2_ndf >= m_pre_track_quality_max ||
       track2.fit_chi2_ndf >= m_pre_track_quality_max))
  {
    return false;
  }

  if (m_pre_track_pt_min > 0.0 &&
      (pt(track1.momentum) < m_pre_track_pt_min || pt(track2.momentum) < m_pre_track_pt_min))
  {
    return false;
  }

  auto dca1 = track1.vertex_dca;
  if (!track1.has_vertex_dca)
  {
    dca1 = fitted_track_dca_to_vertex(track1, primary_vertex);
  }
  auto dca2 = track2.vertex_dca;
  if (!track2.has_vertex_dca)
  {
    dca2 = fitted_track_dca_to_vertex(track2, primary_vertex);
  }
  if (!std::isfinite(dca1.first) || !std::isfinite(dca1.second) ||
      !std::isfinite(dca2.first) || !std::isfinite(dca2.second))
  {
    return false;
  }

  if (m_pre_track_dca_xy_min >= 0.0 &&
      (dca1.first < m_pre_track_dca_xy_min || dca2.first < m_pre_track_dca_xy_min))
  {
    return false;
  }
  if (m_pre_track_dca_z_min >= 0.0 &&
      (dca1.second < m_pre_track_dca_z_min || dca2.second < m_pre_track_dca_z_min))
  {
    return false;
  }
  if (m_pre_track_dca_xy_max >= 0.0 &&
      (dca1.first > m_pre_track_dca_xy_max || dca2.first > m_pre_track_dca_xy_max))
  {
    return false;
  }
  if (m_pre_track_dca_z_max >= 0.0 &&
      (dca1.second > m_pre_track_dca_z_max || dca2.second > m_pre_track_dca_z_max))
  {
    return false;
  }

  if (m_pre_pair_dca_max < 0.0 && m_pre_lproj_min < 0.0 && m_pre_cos_theta_min < -1.0)
  {
    return true;
  }

  LinePca rough_pca;
  if (!line_line_pca(track1.position, track1.momentum, track2.position, track2.momentum, rough_pca, true))
  {
    return false;
  }

  if (m_pre_pair_dca_max >= 0.0 && rough_pca.dca > m_pre_pair_dca_max)
  {
    return false;
  }

  const Vec3 rough_vertex = scale(add(rough_pca.pca1, rough_pca.pca2), 0.5);
  const Vec3 flight = subtract(rough_vertex, primary_vertex);
  const Vec3 total_mom = add(track1.momentum, track2.momentum);
  const double lproj = norm(flight);
  const double cos_theta = vector_cosine(flight, total_mom);

  if (m_pre_lproj_min >= 0.0 && (!std::isfinite(lproj) || lproj < m_pre_lproj_min))
  {
    return false;
  }
  if (m_pre_cos_theta_min >= -1.0 &&
      (!std::isfinite(cos_theta) || cos_theta < m_pre_cos_theta_min))
  {
    return false;
  }

  return true;
}

bool TpcV0CandidateTree::passes_pair_selection(const Vec3 &pca1, const Vec3 &pca2,
                                               const Vec3 &pair_vertex,
                                               const Vec3 &primary_vertex,
                                               const double pair_dca,
                                               const double cos_theta,
                                               const double alpha) const
{
  if (m_pair_pca_z_max >= 0.0 &&
      (!std::isfinite(pair_vertex.z) || std::abs(pair_vertex.z) >= m_pair_pca_z_max))
  {
    return false;
  }

  const double pca_dz = std::abs(pca1.z - pca2.z);
  if (m_pair_pca_dz_max >= 0.0 &&
      (!std::isfinite(pca_dz) || pca_dz >= m_pair_pca_dz_max))
  {
    return false;
  }

  const double dx = pair_vertex.x - primary_vertex.x;
  const double dy = pair_vertex.y - primary_vertex.y;
  const double decay_radius = std::hypot(dx, dy);
  if (m_pair_decay_radius_min >= 0.0 &&
      (!std::isfinite(decay_radius) || decay_radius <= m_pair_decay_radius_min))
  {
    return false;
  }

  if (m_pair_alpha_abs_max >= 0.0 &&
      (!std::isfinite(alpha) || std::abs(alpha) >= m_pair_alpha_abs_max))
  {
    return false;
  }

  if (m_pair_dca_max >= 0.0 &&
      (!std::isfinite(pair_dca) || pair_dca >= m_pair_dca_max))
  {
    return false;
  }

  if (m_pair_dira_min >= -1.0 &&
      (!std::isfinite(cos_theta) || cos_theta <= m_pair_dira_min))
  {
    return false;
  }

  return true;
}
