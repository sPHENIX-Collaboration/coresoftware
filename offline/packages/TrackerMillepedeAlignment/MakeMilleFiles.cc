#include "MakeMilleFiles.h"

#include "Mille.h"

/// Tracking includes

#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TpcDefs.h>  // for side
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <trackreco/ActsPropagator.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <Acts/Definitions/Algebra.hpp>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TF1.h>

#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
  /// square
  template <class T>
  constexpr T square(const T& x)
  {
    return x * x;
  }

  inline bool allFinite(const float* values, int n)
  {
    for (int i = 0; i < n; ++i)
    {
      if (!std::isfinite(values[i]))
      {
        return false;
      }
    }
    return true;
  }
}  // namespace

//____________________________________________________________________________..
MakeMilleFiles::MakeMilleFiles(const std::string& name)
  : SubsysReco(name)
  , _mille(nullptr)
{
}

//____________________________________________________________________________..
int MakeMilleFiles::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  // Instantiate Mille and open output data file
  //  _mille = new Mille(data_outfilename.c_str(), false);   // write text in data files, rather than binary, for debugging only
  _mille = new Mille(data_outfilename.c_str(), _binary);

  // Write the steering file here, and add the data file path to it
  std::ofstream steering_file(steering_outfilename);
  steering_file << data_outfilename << std::endl;
  steering_file << m_constraintFileName << std::endl;
  steering_file.close();

  m_constraintFile.open(m_constraintFileName);
  if (m_useEventVertex)
  {
    for (int i : AlignmentDefs::glbl_vtx_label)
    {
      m_constraintFile << " Constraint   0.0" << std::endl;
      m_constraintFile << "       " << i << "   1" << std::endl;
    }
  }
  // print grouping setup to log file:
  std::cout << "MakeMilleFiles::InitRun: Surface groupings are mvtx " << mvtx_group << " intt " << intt_group << " tpc " << tpc_group << " mms " << mms_group << std::endl;

  // Make TTree for cross check
  if (!m_tfile_name.empty())
  {
    m_file = TFile::Open(m_tfile_name.c_str(), "RECREATE");
    m_trackCutFlow = new TH1D(
        "h_trackCutFlow", "Cumulative MakeMilleFiles track cut flow;selection stage;tracks",
        10, 0.5, 10.5);
    m_trackCutFlow->GetXaxis()->SetBinLabel(1, "State-map entries");
    m_trackCutFlow->GetXaxis()->SetBinLabel(2, "Valid track");
    m_trackCutFlow->GetXaxis()->SetBinLabel(3, "p_{T}");
    m_trackCutFlow->GetXaxis()->SetBinLabel(4, "N_{MVTX}");
    m_trackCutFlow->GetXaxis()->SetBinLabel(5, "N_{INTT}");
    m_trackCutFlow->GetXaxis()->SetBinLabel(6, "N_{vertex tracks}");
    m_trackCutFlow->GetXaxis()->SetBinLabel(7, "Max N_{vertex tracks}(p_{T})");
    m_trackCutFlow->GetXaxis()->SetBinLabel(8, "Vertex association (required)");
    m_trackCutFlow->GetXaxis()->SetBinLabel(9, "Vertex window (required)");
    m_trackCutFlow->GetXaxis()->SetBinLabel(10, "Vertex constraint (required)");
    m_trackCutFlow->SetDirectory(m_file);

    h_track_pt_after_pruning = new TH1D(
        "h_track_pt_after_pruning", "Track p_{T} received by MakeMilleFiles;p_{T} [GeV/c];tracks",
        200, 0, 10);
    h_track_pt_after_pruning->SetDirectory(m_file);

    h_track_pt_after_makemille_pt_cut = new TH1D(
        "h_track_pt_after_makemille_pt_cut", "Track p_{T} after MakeMilleFiles pT cut;p_{T} [GeV/c];tracks",
        200, 0, 10);
    h_track_pt_after_makemille_pt_cut->SetDirectory(m_file);

    h_vertex_ntracks = new TH1D(
        "h_vertex_ntracks", "Vertex track count used by the MakeMilleFiles cut (-2 = no vertex, -1 = crossing mismatch);vertex_ntracks;tracks",
        103, -2.5, 100.5);
    h_vertex_ntracks->SetDirectory(m_file);

    h_vertex_ntracks_raw = new TH1D(
        "h_vertex_ntracks_raw", "Track count of the assigned vertex, ignoring crossing match;vertex tracks;tracks",
        101, -0.5, 100.5);
    h_vertex_ntracks_raw->SetDirectory(m_file);

    // The last bin of h_absent_vertex_id and the overflow of this histogram hold
    // the unassigned (maximum unsigned) vertex ID.
    h_absent_vertex_id = new TH1D(
        "h_absent_vertex_id", "Requested vertex ID of tracks with no vertex in the map (overflow = unassigned ID);vertex ID;tracks",
        101, -0.5, 100.5);
    h_absent_vertex_id->SetDirectory(m_file);

    h_track_crossing_with_vertex = new TH1D(
        "h_track_crossing_with_vertex", "Track crossing, tracks whose vertex exists;track crossing;tracks",
        1001, -500.5, 500.5);
    h_track_crossing_with_vertex->SetDirectory(m_file);

    h_track_crossing_without_vertex = new TH1D(
        "h_track_crossing_without_vertex", "Track crossing, tracks whose vertex is absent;track crossing;tracks",
        1001, -500.5, 500.5);
    h_track_crossing_without_vertex->SetDirectory(m_file);

    h_track_vertex_crossing_diff = new TH1D(
        "h_track_vertex_crossing_diff", "Track crossing - vertex crossing (tracks with a vertex);track crossing - vertex crossing;tracks",
        101, -50.5, 50.5);
    h_track_vertex_crossing_diff->SetDirectory(m_file);

    h_frac_no_vertex_vs_ntracks = new TH2D(
        "h_frac_no_vertex_vs_ntracks",
        "Per crossing: fraction of tracks without vertex;tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);fraction without vertex",
        100, 0.5, 100.5, 21, -0.025, 1.025);
    h_frac_no_vertex_vs_ntracks->SetDirectory(m_file);

    p_frac_no_vertex_vs_ntracks = new TProfile(
        "p_frac_no_vertex_vs_ntracks",
        "Per crossing: mean fraction of tracks without vertex;tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);mean fraction without vertex",
        100, 0.5, 100.5);
    p_frac_no_vertex_vs_ntracks->SetDirectory(m_file);

    // Tracks with unknown crossing (SHRT_MAX) are pooled per event.
    h_ntracks_nocrossing = new TH1D(
        "h_ntracks_nocrossing",
        "Tracks with unknown crossing per event (after p_{T}/N_{MVTX}/N_{INTT} cuts);tracks with unknown crossing;events",
        101, -0.5, 100.5);
    h_ntracks_nocrossing->SetDirectory(m_file);

    h_frac_no_vertex_vs_ntracks_nocrossing = new TH2D(
        "h_frac_no_vertex_vs_ntracks_nocrossing",
        "Per event, unknown-crossing tracks: fraction without vertex;tracks with unknown crossing in event;fraction without vertex",
        100, 0.5, 100.5, 21, -0.025, 1.025);
    h_frac_no_vertex_vs_ntracks_nocrossing->SetDirectory(m_file);

    p_frac_no_vertex_vs_ntracks_nocrossing = new TProfile(
        "p_frac_no_vertex_vs_ntracks_nocrossing",
        "Per event, unknown-crossing tracks: mean fraction without vertex;tracks with unknown crossing in event;mean fraction without vertex",
        100, 0.5, 100.5);
    p_frac_no_vertex_vs_ntracks_nocrossing->SetDirectory(m_file);

    // vertex_ntracks of the assigned vertex, crossing-matched tracks only
    h_vertex_ntracks_vs_ntracks_crossing = new TH2D(
        "h_vertex_ntracks_vs_ntracks_crossing",
        "Per track: assigned vertex track count vs tracks in crossing (crossing-matched tracks);tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);vertex_ntracks (crossing matched)",
        100, 0.5, 100.5, 100, 0.5, 100.5);
    h_vertex_ntracks_vs_ntracks_crossing->SetDirectory(m_file);

    // Same tracks as h_vertex_ntracks_vs_ntracks_crossing (crossing matched)
    p_vertex_ntracks_vs_ntracks_crossing = new TProfile(
        "p_vertex_ntracks_vs_ntracks_crossing",
        "Per track: mean assigned vertex track count vs tracks in crossing (crossing-matched tracks);tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);mean vertex_ntracks (crossing matched)",
        100, 0.5, 100.5);
    p_vertex_ntracks_vs_ntracks_crossing->SetDirectory(m_file);

    // Split of the crossing-matched tracks: used to build the vertex
    // (in vertex->begin_tracks()...end_tracks()) or attached afterwards by
    // PHSimpleVertexFinder's nearest-z step (not in the vertex track list)
    h_vertex_ntracks_vs_ntracks_crossing_used = new TH2D(
        "h_vertex_ntracks_vs_ntracks_crossing_used",
        "Per track: assigned vertex track count vs tracks in crossing (tracks used to build the vertex);tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);vertex_ntracks (crossing matched)",
        100, 0.5, 100.5, 100, 0.5, 100.5);
    h_vertex_ntracks_vs_ntracks_crossing_used->SetDirectory(m_file);

    p_vertex_ntracks_vs_ntracks_crossing_used = new TProfile(
        "p_vertex_ntracks_vs_ntracks_crossing_used",
        "Per track: mean assigned vertex track count vs tracks in crossing (tracks used to build the vertex);tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);mean vertex_ntracks (crossing matched)",
        100, 0.5, 100.5);
    p_vertex_ntracks_vs_ntracks_crossing_used->SetDirectory(m_file);

    h_vertex_ntracks_vs_ntracks_crossing_attached = new TH2D(
        "h_vertex_ntracks_vs_ntracks_crossing_attached",
        "Per track: assigned vertex track count vs tracks in crossing (tracks attached after vertexing);tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);vertex_ntracks (crossing matched)",
        100, 0.5, 100.5, 100, 0.5, 100.5);
    h_vertex_ntracks_vs_ntracks_crossing_attached->SetDirectory(m_file);

    p_vertex_ntracks_vs_ntracks_crossing_attached = new TProfile(
        "p_vertex_ntracks_vs_ntracks_crossing_attached",
        "Per track: mean assigned vertex track count vs tracks in crossing (tracks attached after vertexing);tracks in crossing (after p_{T}/N_{MVTX}/N_{INTT} cuts);mean vertex_ntracks (crossing matched)",
        100, 0.5, 100.5);
    p_vertex_ntracks_vs_ntracks_crossing_attached->SetDirectory(m_file);

    m_ntuple = new TNtuple (
      "ntp", "ntp",
      "layer:stave:chip:residualX:residualY:clusphi:xglob:yglob:zglob:"
      "errX:errY:"
      "dXdR:dXdZ0:dXdphi:dXdtheta:dXdqoverp:dXdt:"
      "dXdalpha:dXdbeta:dXdgamma:dXdx:dXdy:dXdz:"
      "dYdR:dYdZ0:dYdphi:dYdtheta:dYdqoverp:dYdt:"
      "dYdalpha:dYdbeta:dYdgamma:dYdx:dYdy:dYdz"
    );
    m_ntuple->SetDirectory(m_file);

    track_ntp = new TNtuple(
        "track_ntp", "MakeMilleFiles track ntuple",
        "track_id:quality:residual_vertexX:residual_vertexY:"
        "residualxsigma:residualysigma:"
        "dXdR:dXdZ0:dXdphi:dXdtheta:dXdqoverp:dXdt:"
        "dYdR:dYdZ0:dYdphi:dYdtheta:dYdqoverp:dYdt:"
        "dXdx:dXdy:dXdz:dYdx:dYdy:dYdz:"
        "track_xvtx:track_yvtx:track_zvtx:"
        "event_xvtx:event_yvtx:event_zvtx:track_phi:track_eta:track_p:track_pt:"
        "n_tracks_vtx");
    track_ntp->SetDirectory(m_file);
  }

  return ret;
}

//____________________________________________________________________________..
int MakeMilleFiles::process_event(PHCompositeNode* /*topNode*/)
{
  // Outline:
  //
  // loop over track alignment states
  //   Make any track cuts here to skip undesirable tracks (maybe low pT?)
  //   loop over track states+measurements for each track
  //      for each measurement, performed in trackreco/ActsAlignmentStates.cc
  //         Get measurement value and error
  //         Calculate derivatives and residuals from Acts jacobians
  //         These are stored in a map and unpacked for mille
  //   Call _mille->mille() with arguments obtained from previous iteration:
  //     local pars
  //     array of local derivatives
  //     global pars
  //     array of global derivatives
  //     array of integer global par labels
  //     residual value (float) z = measurement - track state
  //     sigma of measurement
  //   After processing all measurements for this track, call _mille->end() to add buffer to file and reset buffer
  // After all tracks are processed, file is closed when Mille destructor is called
  // Note: all units are in the Acts units of mm and GeV to avoid converting matrices

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;
    std::cout << "state map size " << _state_map->size() << std::endl;
  }

  ActsPropagator propagator(_tGeometry);

  unsigned int n_track_map_matches = 0;
  unsigned int n_track_cuts_pass = 0;
  unsigned int n_matching_vertices = 0;
  unsigned int n_vertex_position_matches = 0;
  unsigned int n_propagation_successes = 0;
  unsigned int n_finite_vertex_residuals = 0;
  unsigned int n_track_ntp_fills = 0;
  unsigned int n_missing_track_map_entries = 0;
  unsigned int n_null_tracks = 0;
  unsigned int n_empty_state_vectors = 0;
  unsigned int n_rejected_pt = 0;
  unsigned int n_rejected_nmvtx = 0;
  unsigned int n_rejected_nintt = 0;
  unsigned int n_rejected_vertex_multiplicity = 0;
  unsigned int n_rejected_vertex_high_pt_multiplicity = 0;
  unsigned int n_missing_matching_vertices = 0;
  unsigned int n_track_vertex_absent = 0;
  unsigned int n_track_vertex_crossing_mismatch = 0;
  unsigned int n_absent_vertex_id_unassigned = 0;
  unsigned int n_crossing_diff_underflow = 0;
  unsigned int n_crossing_diff_overflow = 0;
  unsigned int n_rejected_vertex_position = 0;
  unsigned int n_insufficient_track_states = 0;
  unsigned int n_propagation_failures = 0;
  unsigned int n_nonfinite_vertex_residuals = 0;
  unsigned int n_unavailable_track_ntuple = 0;
  std::unordered_map<unsigned int, unsigned int> high_pt_tracks_per_vertex;
  bool printed_vertex_map = false;
  // crossing -> (tracks passing pT/nMVTX/nINTT cuts, of which without a vertex)
  std::map<short int, std::pair<unsigned int, unsigned int>> tracks_per_crossing;
  // Same counts for tracks with unknown crossing (SHRT_MAX), pooled per event
  unsigned int n_nocrossing_tracks = 0;
  unsigned int n_nocrossing_without_vertex = 0;
  // crossing -> vertex->size_tracks() of crossing-matched tracks only
  std::map<short int, std::vector<int>> matched_vertex_ntracks_per_crossing;
  // Same, split into tracks used to build the vertex and tracks attached afterwards
  std::map<short int, std::vector<int>> used_vertex_ntracks_per_crossing;
  std::map<short int, std::vector<int>> attached_vertex_ntracks_per_crossing;

  for (auto [key, statevec] : *_state_map)
  {
    if (m_trackCutFlow)
    {
      m_trackCutFlow->Fill(1.0);
    }

    // Check if track was removed from cleaner
    auto iter = _track_map->find(key);
    if (iter == _track_map->end())
    {
      ++n_missing_track_map_entries;
      if (m_trackNtpDiagnostics)
      {
        std::cout << "MakeMilleFiles track_ntp reject: alignment-state key "
                  << key << " is absent from " << m_track_map_name
                  << std::endl;
      }
      continue;
    }

    SvtxTrack* track = iter->second;
    if (!track)
    {
      ++n_null_tracks;
      if (m_trackNtpDiagnostics)
      {
        std::cout << "MakeMilleFiles track_ntp reject: null track for key "
                  << key << std::endl;
      }
      continue;
    }
    ++n_track_map_matches;
    if (m_trackCutFlow)
    {
      m_trackCutFlow->Fill(2.0);
    }
    if (h_track_pt_after_pruning)
    {
      h_track_pt_after_pruning->Fill(track->get_pt());
    }

    if (statevec.empty())
    {
      ++n_empty_state_vectors;
      if (m_trackNtpDiagnostics)
      {
        std::cout << "MakeMilleFiles track_ntp diagnostic: track "
                  << track->get_id()
                  << " has an empty alignment-state vector" << std::endl;
      }
    }

    if (Verbosity() > 0)
    {
      std::cout << std::endl
                << __LINE__ << ": Processing track itrack: " << key << ": nhits: " << track->size_cluster_keys()
                << ": Total tracks: " << _track_map->size() << ": phi: " << track->get_phi() << std::endl;
    }

    // Apply the track-quality cuts before adding any cluster or vertex residuals.
    unsigned int nmvtx = 0;
    unsigned int nintt = 0;
    const auto count_silicon_clusters = [&nmvtx, &nintt](const auto* cluster_source)
    {
      for (auto cluster_iter = cluster_source->begin_cluster_keys();
           cluster_iter != cluster_source->end_cluster_keys(); ++cluster_iter)
      {
        switch (TrkrDefs::getTrkrId(*cluster_iter))
        {
        case TrkrDefs::mvtxId:
          ++nmvtx;
          break;

        case TrkrDefs::inttId:
          ++nintt;
          break;

        default:
          break;
        }
      }
    };

    if (track->begin_cluster_keys() != track->end_cluster_keys())
    {
      count_silicon_clusters(track);
    }
    else if (const TrackSeed* silicon_seed = track->get_silicon_seed())
    {
      count_silicon_clusters(silicon_seed);
    }

    size_t vertex_ntracks = 0;
    unsigned int n_vertex_tracks_above_pt = 0;
    const bool applyMaxTracksAboveMinPtCut =
        m_maxTracksAboveMinPtPerVertex !=
        std::numeric_limits<unsigned int>::max();
    const SvtxVertex* vertex = nullptr;
    if (m_useEventVertex && m_vertex_map)
    {
      vertex = m_vertex_map->get(track->get_vertex_id());
      if (!vertex)
      {
        ++n_track_vertex_absent;
        if (track->get_vertex_id() == std::numeric_limits<unsigned int>::max())
        {
          ++n_absent_vertex_id_unassigned;
        }
        if (h_absent_vertex_id)
        {
          h_absent_vertex_id->Fill(static_cast<double>(track->get_vertex_id()));
        }
        if (h_track_crossing_without_vertex)
        {
          h_track_crossing_without_vertex->Fill(track->get_crossing());
        }
      }
      else
      {
        const int crossing_diff =
            static_cast<int>(track->get_crossing()) -
            static_cast<int>(vertex->get_beam_crossing());
        if (crossing_diff != 0)
        {
          ++n_track_vertex_crossing_mismatch;
        }
        if (crossing_diff < -50)
        {
          ++n_crossing_diff_underflow;
        }
        else if (crossing_diff > 50)
        {
          ++n_crossing_diff_overflow;
        }
        if (h_vertex_ntracks_raw)
        {
          h_vertex_ntracks_raw->Fill(vertex->size_tracks());
        }
        if (h_track_crossing_with_vertex)
        {
          h_track_crossing_with_vertex->Fill(track->get_crossing());
        }
        if (h_track_vertex_crossing_diff)
        {
          h_track_vertex_crossing_diff->Fill(crossing_diff);
        }
      }
      if (vertex && vertex->get_beam_crossing() == track->get_crossing())
      {
        vertex_ntracks = vertex->size_tracks();
        if (applyMaxTracksAboveMinPtCut)
        {
          const unsigned int vertex_id = track->get_vertex_id();
          const auto cached_count = high_pt_tracks_per_vertex.find(vertex_id);
          if (cached_count != high_pt_tracks_per_vertex.end())
          {
            n_vertex_tracks_above_pt = cached_count->second;
          }
          else
          {
            for (auto vertex_track_iter = vertex->begin_tracks();
                 vertex_track_iter != vertex->end_tracks(); ++vertex_track_iter)
            {
              const auto track_iter = _track_map->find(*vertex_track_iter);
              if (track_iter != _track_map->end() &&
                  track_iter->second &&
                  track_iter->second->get_pt() >= m_minPt)
              {
                ++n_vertex_tracks_above_pt;
              }
            }
            high_pt_tracks_per_vertex.emplace(vertex_id,
                                              n_vertex_tracks_above_pt);
          }
        }
      }
      if (h_vertex_ntracks)
      {
        // -2: no vertex in the map, -1: vertex crossing differs from the track's
        double vertex_ntracks_fill = vertex_ntracks;
        if (!vertex)
        {
          vertex_ntracks_fill = -2;
        }
        else if (vertex->get_beam_crossing() != track->get_crossing())
        {
          vertex_ntracks_fill = -1;
        }
        h_vertex_ntracks->Fill(vertex_ntracks_fill);
      }
    }

    if (m_useEventVertex && Verbosity() > 0)
    {
      if (vertex && vertex->get_beam_crossing() == track->get_crossing())
      {
        std::cout << "Track " << track->get_id()
                  << " matched vertex " << track->get_vertex_id()
                  << " at crossing " << track->get_crossing()
                  << ", position (cm) = (" << vertex->get_x() << ", "
                  << vertex->get_y() << ", " << vertex->get_z() << ")"
                  << std::endl;
      }
      else
      {
        std::cout << "Track " << track->get_id()
                  << " has no matching vertex: requested vertex "
                  << track->get_vertex_id() << " at crossing "
                  << track->get_crossing();
        if (vertex)
        {
          std::cout << "; assigned vertex crossing "
                    << vertex->get_beam_crossing()
                    << ", position (cm) = (" << vertex->get_x() << ", "
                    << vertex->get_y() << ", " << vertex->get_z() << ")";
        }
        else
        {
          std::cout << "; requested vertex is absent";
        }
        std::cout << std::endl;

        if (!printed_vertex_map && m_vertex_map)
        {
          std::cout << "Available vertices:" << std::endl;
          for (const auto& [vertex_id, available_vertex] : *m_vertex_map)
          {
            if (available_vertex)
            {
              std::cout << "  vertex " << vertex_id
                        << ", crossing "
                        << available_vertex->get_beam_crossing()
                        << ", position (cm) = ("
                        << available_vertex->get_x() << ", "
                        << available_vertex->get_y() << ", "
                        << available_vertex->get_z() << ")"
                        << ", ntracks=" << available_vertex->size_tracks()
                        << std::endl;
            }
          }
          printed_vertex_map = true;
        }
      }
    }

    const bool hasTooManyTracksAboveMinPt =
        applyMaxTracksAboveMinPtCut &&
        n_vertex_tracks_above_pt > m_maxTracksAboveMinPtPerVertex;
    const bool failsPt = track->get_pt() < m_minPt;
    if (h_track_pt_after_makemille_pt_cut && !failsPt)
    {
      h_track_pt_after_makemille_pt_cut->Fill(track->get_pt());
    }
    const bool failsNmvtx = nmvtx < m_minNmvtx;
    const bool failsNintt = nintt < m_minNintt;
    const bool failsVertexMultiplicity =
        vertex_ntracks < m_minVertexTracks;

    // Crossing-by-crossing vertex association, before any vertex requirement.
    // "Without vertex" uses the same test as getEventVertex() (cut-flow bin 8).
    if (m_useEventVertex && m_vertex_map && !failsPt && !failsNmvtx &&
        !failsNintt)
    {
      const bool hasVertex =
          vertex && vertex->get_beam_crossing() == track->get_crossing();
      if (track->get_crossing() != SHRT_MAX)
      {
        auto& [n_tracks, n_without_vertex] =
            tracks_per_crossing[track->get_crossing()];
        ++n_tracks;
        n_without_vertex += !hasVertex;
        // For h_/p_vertex_ntracks_vs_ntracks_crossing: only tracks whose vertex
        // has the same crossing (same test as cut-flow bin 8).
        if (hasVertex)
        {
          matched_vertex_ntracks_per_crossing[track->get_crossing()].push_back(
              static_cast<int>(vertex->size_tracks()));
          // key is the track-map key, which PHSimpleVertexFinder stores in the vertex
          auto& role_list = (vertex->find_track(key) != vertex->end_tracks())
                                ? used_vertex_ntracks_per_crossing[track->get_crossing()]
                                : attached_vertex_ntracks_per_crossing[track->get_crossing()];
          role_list.push_back(static_cast<int>(vertex->size_tracks()));
        }
      }
      else
      {
        // Unknown crossing: PHSimpleVertexFinder vertexes these tracks as their
        // own group, so a matching vertex also has crossing SHRT_MAX.
        ++n_nocrossing_tracks;
        n_nocrossing_without_vertex += !hasVertex;
      }
    }
    if (m_trackCutFlow && !failsPt)
    {
      m_trackCutFlow->Fill(3.0);
    }
    if (m_trackCutFlow && !failsPt && !failsNmvtx)
    {
      m_trackCutFlow->Fill(4.0);
    }
    if (m_trackCutFlow && !failsPt && !failsNmvtx && !failsNintt)
    {
      m_trackCutFlow->Fill(5.0);
    }
    if (m_trackCutFlow && !failsPt && !failsNmvtx && !failsNintt &&
        !failsVertexMultiplicity)
    {
      m_trackCutFlow->Fill(6.0);
    }
    if (m_trackCutFlow && !failsPt && !failsNmvtx && !failsNintt &&
        !failsVertexMultiplicity && !hasTooManyTracksAboveMinPt)
    {
      m_trackCutFlow->Fill(7.0);
    }
    if (failsPt || failsNmvtx || failsNintt || failsVertexMultiplicity ||
        hasTooManyTracksAboveMinPt)
    {
      n_rejected_pt += failsPt;
      n_rejected_nmvtx += failsNmvtx;
      n_rejected_nintt += failsNintt;
      n_rejected_vertex_multiplicity += failsVertexMultiplicity;
      n_rejected_vertex_high_pt_multiplicity +=
          hasTooManyTracksAboveMinPt;
      if (Verbosity() > 0)
      {
        std::cout << "Skipping track with pT=" << track->get_pt()
                  << ", nMVTX=" << nmvtx << ", nINTT=" << nintt
                  << ", vertex_ntracks=" << vertex_ntracks;
        if (applyMaxTracksAboveMinPtCut)
        {
          std::cout << ", vertex tracks with pT >= " << m_minPt
                    << "=" << n_vertex_tracks_above_pt
                    << ", maximum=" << m_maxTracksAboveMinPtPerVertex;
        }
        std::cout << " (minimums: " << m_minPt << ", "
                  << m_minNmvtx << ", " << m_minNintt << ", "
                  << m_minVertexTracks << ")"
                  << std::endl;
      }
      if (m_trackNtpDiagnostics)
      {
        std::cout << "MakeMilleFiles track_ntp reject: track "
                  << track->get_id() << " failed"
                  << (failsPt ? " pT" : "")
                  << (failsNmvtx ? " nMVTX" : "")
                  << (failsNintt ? " nINTT" : "")
                  << (failsVertexMultiplicity ? " vertex_ntracks: " : "")
                  <<vertex_ntracks<<" "
                  << (hasTooManyTracksAboveMinPt
                          ? " vertex_high_pt_multiplicity"
                          : "")
                  << std::endl;
      }
      continue;
    }
    ++n_track_cuts_pass;

    // With the event vertex enabled, a track must have a matching vertex
    // inside the configured x/y/z window before any of its residuals are used.
    Acts::Vector3 eventVertex = Acts::Vector3::Zero();
    const SvtxVertex* svtxVertex = nullptr;
    bool isVertexPositionMatch = false;
    if (m_useEventVertex)
    {
      eventVertex = getEventVertex(track, svtxVertex);
      if (!svtxVertex || !std::isfinite(eventVertex.z()))
      {
        ++n_missing_matching_vertices;
        if (m_trackNtpDiagnostics)
        {
          std::cout << "MakeMilleFiles track_ntp reject: track "
                    << track->get_id() << " has no vertex matching ID "
                    << track->get_vertex_id() << " and crossing "
                    << track->get_crossing() << std::endl;
        }
        continue;
      }
      ++n_matching_vertices;
      if (m_trackCutFlow)
      {
        m_trackCutFlow->Fill(8.0);
      }

      isVertexPositionMatch =
          std::abs(track->get_x() - eventVertex.x()) < m_vertexMatchCutCm.x() &&
          std::abs(track->get_y() - eventVertex.y()) < m_vertexMatchCutCm.y() &&
          std::abs(track->get_z() - eventVertex.z()) < m_vertexMatchCutCm.z();
      if (!isVertexPositionMatch)
      {
        ++n_rejected_vertex_position;
        if (m_trackNtpDiagnostics)
        {
          std::cout << "MakeMilleFiles track_ntp reject: track "
                    << track->get_id() << " is outside the vertex window; "
                    << "track (cm) = (" << track->get_x() << ", "
                    << track->get_y() << ", " << track->get_z()
                    << "), vertex (cm) = (" << eventVertex.x() << ", "
                    << eventVertex.y() << ", " << eventVertex.z() << ")"
                    << std::endl;
        }
        continue;
      }
      ++n_vertex_position_matches;
      if (m_trackCutFlow)
      {
        m_trackCutFlow->Fill(9.0);
      }
    }

    // With the event vertex, the cluster residuals are written only after the
    // vertex constraint is verified to be usable (see below).
    if (!m_useEventVertex)
    {
      addTrackToMilleFile(statevec);
    }
    else
    {
      //! Retain the reconstructed transverse vertex coordinates.
      // eventVertex(0) = 0;
      // eventVertex(1) = 0;

      auto dcapair = TrackAnalysisUtils::get_dca(track, eventVertex);
      Acts::Vector2 vtx_residual(-dcapair.first.first, -dcapair.second.first);
      // Convert the DCA residual from cm to mm; m_vtxSigma is specified in mm.
      vtx_residual *= Acts::UnitConstants::cm;

      float lclvtx_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NLOC];
      Acts::BoundSquareMatrix perigeeCov = Acts::BoundSquareMatrix::Zero();
      float lclvtx_derivative_psuedo[SvtxAlignmentState::NLOC][SvtxAlignmentState::NLOC]{};
      float localmeaserrpsuedo[SvtxAlignmentState::NLOC]{};
      bool success = false;
      bool success2 = false;
      if (track->size_states() < 2)
      {
        std::cout<<"track->size_states(): "<<track->size_states()<<std::endl;
        ++n_insufficient_track_states;
        if (m_trackNtpDiagnostics)
        {
          std::cout << "MakeMilleFiles track_ntp reject: track "
                    << track->get_id() << " has " << track->size_states()
                    << " states; need at least two to propagate from a "
                       "silicon measurement"
                    << std::endl;
        }
      }
      else
      {
        success = getLocalVtxDerivativesXY(track, propagator,
                                           eventVertex, lclvtx_derivative,
                                           &perigeeCov);
        // Pseudo-measurements need the track covariance on the vertex
        // perigee surface, so only build them after a successful propagation.
        if (success)
        {
          success2 = getLocalVtxDerivativesXYpsuedo(track, svtxVertex, perigeeCov,
                                                    lclvtx_derivative_psuedo,
                                                    localmeaserrpsuedo);
          std::pair<float,float> rotatedVtxCovDiag = getRotatedVtxCovarianceDiagonalPerigee(track, svtxVertex);
          set_vtx_sigma(rotatedVtxCovDiag.first, rotatedVtxCovDiag.second);
        }
      }

      // The global derivs dimensions are [alpha/beta/gamma](x/y/z)
      float glblvtx_derivative[SvtxAlignmentState::NRES][3];
      getGlobalVtxDerivativesXY(track, eventVertex, glblvtx_derivative);

      if (Verbosity() > 2)
      {
        std::cout << "vertex info for trakc " << track->get_id() << " with charge " << track->get_charge() << std::endl;
        std::cout << "vertex is " << eventVertex.transpose() << std::endl;
        std::cout << "vertex residuals " << vtx_residual.transpose()
                  << std::endl;
        std::cout << "global vtx derivatives " << std::endl;
        for (auto& i : glblvtx_derivative)
        {
          for (float j : i)
          {
            std::cout << j << ", ";
          }
          std::cout << std::endl;
        }
      }
      // track_ntp records only tracks with successful vertex propagation and
      // finite vertex residuals and derivatives in both components.
      if (success)
      {
        ++n_propagation_successes;
        // The vertex constraint is used only if both residual components and
        // all of their local and global derivatives are finite.
        bool hasFiniteVertexResidual = true;
        for (int i = 0; i < 2; i++)
        {
          hasFiniteVertexResidual =
              hasFiniteVertexResidual && std::isfinite(vtx_residual(i)) &&
              allFinite(lclvtx_derivative[i], SvtxAlignmentState::NLOC) &&
              allFinite(glblvtx_derivative[i], AlignmentDefs::NGLVTX);
        }
        if (!hasFiniteVertexResidual)
        {
          ++n_nonfinite_vertex_residuals;
          if (m_trackNtpDiagnostics)
          {
            std::cout << "MakeMilleFiles track_ntp reject: track "
                      << track->get_id()
                      << " has non-finite vertex residuals or derivatives"
                      << std::endl;
          }
          continue;
        }
        ++n_finite_vertex_residuals;
        if (m_trackCutFlow)
        {
          m_trackCutFlow->Fill(10.0);
        }
        addTrackToMilleFile(statevec);

        if (!track_ntp)
        {
          ++n_unavailable_track_ntuple;
          if (m_trackNtpDiagnostics)
          {
            std::cout << "MakeMilleFiles track_ntp reject: track_ntp tree "
                      << "is unavailable for track " << track->get_id()
                      << std::endl;
          }
        }
        else
        {
          float track_ntp_data[] = {
              static_cast<float>(track->get_id()),
              static_cast<float>(track->get_quality()),
              static_cast<float>(vtx_residual(0)),
              static_cast<float>(vtx_residual(1)),
              static_cast<float>(m_vtxSigma(0) * m_vtxSigmaInflation),
              static_cast<float>(m_vtxSigma(1) * m_vtxSigmaInflation),
              lclvtx_derivative[0][0], lclvtx_derivative[0][1],
              lclvtx_derivative[0][2], lclvtx_derivative[0][3],
              lclvtx_derivative[0][4], lclvtx_derivative[0][5],
              lclvtx_derivative[1][0], lclvtx_derivative[1][1],
              lclvtx_derivative[1][2], lclvtx_derivative[1][3],
              lclvtx_derivative[1][4], lclvtx_derivative[1][5],
              glblvtx_derivative[0][0], glblvtx_derivative[0][1],
              glblvtx_derivative[0][2],
              glblvtx_derivative[1][0], glblvtx_derivative[1][1],
              glblvtx_derivative[1][2],
              // SvtxTrack and SvtxVertex positions are in cm; store mm to
              // match the vertex residual and m_vtxSigma conventions.
              static_cast<float>(track->get_x() * Acts::UnitConstants::cm),
              static_cast<float>(track->get_y() * Acts::UnitConstants::cm),
              static_cast<float>(track->get_z() * Acts::UnitConstants::cm),
              static_cast<float>(eventVertex.x() * Acts::UnitConstants::cm),
              static_cast<float>(eventVertex.y() * Acts::UnitConstants::cm),
              static_cast<float>(eventVertex.z() * Acts::UnitConstants::cm),
              static_cast<float>(track->get_phi()),
              static_cast<float>(track->get_eta()),
              static_cast<float>(track->get_p()),
              static_cast<float>(track->get_pt()),
              static_cast<float>(vertex_ntracks)};
          track_ntp->Fill(track_ntp_data);
          ++n_track_ntp_fills;
          if (m_trackNtpDiagnostics)
          {
            std::cout << "MakeMilleFiles track_ntp fill: track "
                      << track->get_id() << ", vertex residuals (mm) = ("
                      << vtx_residual(0) << ", " << vtx_residual(1)
                      << ")" << std::endl;
          }
        }

        if (Verbosity() > 2)
        {
          // Sign check: d(d0)/d(loc0) should be +/-1; compare its sign with the
          // get_dca-based vertex residual to verify the two conventions agree.
          std::cout << "vertex sign check for track " << track->get_id()
                    << ": vertex residual (mm) = (" << vtx_residual(0) << ", "
                    << vtx_residual(1) << "), d(d0)/d(loc0) = "
                    << lclvtx_derivative[0][0] << std::endl;
        }
        for (int i = 0; i < 2; i++)
        {
          _mille->mille(SvtxAlignmentState::NLOC, lclvtx_derivative[i],
                        AlignmentDefs::NGLVTX, glblvtx_derivative[i],
                        AlignmentDefs::glbl_vtx_label, vtx_residual(i),
                        m_vtxSigma(i) * m_vtxSigmaInflation);
        }
        // Vertex pseudo-measurements carry no global dependence
        float glbl_derivative_vtx_dummy_empty[AlignmentDefs::NGLVTX]{};
        for (int i = 0; success2 && i < SvtxAlignmentState::NLOC; ++i)
        {
          if (!(localmeaserrpsuedo[i] > 0) || !std::isfinite(localmeaserrpsuedo[i]) ||
              !allFinite(lclvtx_derivative_psuedo[i], SvtxAlignmentState::NLOC))
          {
            continue;
          }
          _mille->mille(SvtxAlignmentState::NLOC, lclvtx_derivative_psuedo[i],
                        AlignmentDefs::NGLVTX, glbl_derivative_vtx_dummy_empty,
                        AlignmentDefs::glbl_vtx_label, 0.0,
                        localmeaserrpsuedo[i]);
        }
      }
      else
      {
        if (track->size_states() >= 2)
        {
          ++n_propagation_failures;
          if (m_trackNtpDiagnostics)
          {
            std::cout << "MakeMilleFiles track_ntp reject: propagation to "
                      << "vertex failed for track " << track->get_id()
                      << " (vertex ID " << track->get_vertex_id() << ")"
                      << std::endl;
          }
        }
        continue;
      }
    }

    //! Finish this track
    _mille->end();
  }

  for (const auto& crossing_counts : tracks_per_crossing)
  {
    const auto [n_tracks, n_without_vertex] = crossing_counts.second;
    const double fraction =
        static_cast<double>(n_without_vertex) / n_tracks;
    if (h_frac_no_vertex_vs_ntracks)
    {
      h_frac_no_vertex_vs_ntracks->Fill(n_tracks, fraction);
    }
    if (p_frac_no_vertex_vs_ntracks)
    {
      p_frac_no_vertex_vs_ntracks->Fill(n_tracks, fraction);
    }
  }

  for (const auto& [crossing, matched_vertex_ntracks_list] : matched_vertex_ntracks_per_crossing)
  {
    const unsigned int n_tracks = tracks_per_crossing[crossing].first;
    for (const int vtx_ntracks : matched_vertex_ntracks_list)
    {
      if (h_vertex_ntracks_vs_ntracks_crossing)
      {
        h_vertex_ntracks_vs_ntracks_crossing->Fill(n_tracks, vtx_ntracks);
      }
      if (p_vertex_ntracks_vs_ntracks_crossing)
      {
        p_vertex_ntracks_vs_ntracks_crossing->Fill(n_tracks, vtx_ntracks);
      }
    }
  }

  const auto fill_vertex_ntracks =
      [&tracks_per_crossing](const std::map<short int, std::vector<int>>& per_crossing,
                             TH2D* h2, TProfile* prof)
  {
    for (const auto& [crossing, vertex_ntracks_list] : per_crossing)
    {
      const unsigned int n_tracks = tracks_per_crossing[crossing].first;
      for (const int vtx_ntracks : vertex_ntracks_list)
      {
        if (h2)
        {
          h2->Fill(n_tracks, vtx_ntracks);
        }
        if (prof)
        {
          prof->Fill(n_tracks, vtx_ntracks);
        }
      }
    }
  };
  fill_vertex_ntracks(used_vertex_ntracks_per_crossing,
                      h_vertex_ntracks_vs_ntracks_crossing_used,
                      p_vertex_ntracks_vs_ntracks_crossing_used);
  fill_vertex_ntracks(attached_vertex_ntracks_per_crossing,
                      h_vertex_ntracks_vs_ntracks_crossing_attached,
                      p_vertex_ntracks_vs_ntracks_crossing_attached);

  if (m_useEventVertex && m_vertex_map)
  {
    if (h_ntracks_nocrossing)
    {
      h_ntracks_nocrossing->Fill(n_nocrossing_tracks);
    }
    if (n_nocrossing_tracks > 0)
    {
      const double fraction =
          static_cast<double>(n_nocrossing_without_vertex) / n_nocrossing_tracks;
      if (h_frac_no_vertex_vs_ntracks_nocrossing)
      {
        h_frac_no_vertex_vs_ntracks_nocrossing->Fill(n_nocrossing_tracks, fraction);
      }
      if (p_frac_no_vertex_vs_ntracks_nocrossing)
      {
        p_frac_no_vertex_vs_ntracks_nocrossing->Fill(n_nocrossing_tracks, fraction);
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "Finished processing mille file " << std::endl;
    std::cout << "MakeMilleFiles::process_event vertex diagnostics: "
              << "state-map tracks=" << _state_map->size()
              << ", track-map matches=" << n_track_map_matches
              << ", track-cuts-pass=" << n_track_cuts_pass
              << ", matching vertices=" << n_matching_vertices
              << ", within 2 mm=" << n_vertex_position_matches
              << ", propagation successes=" << n_propagation_successes
              << ", finite residuals=" << n_finite_vertex_residuals
              << ", track_ntp fills=" << n_track_ntp_fills
              << ", missing track-map entries=" << n_missing_track_map_entries
              << ", null tracks=" << n_null_tracks
              << ", empty state vectors=" << n_empty_state_vectors
              << ", rejected pT=" << n_rejected_pt
              << ", rejected nMVTX=" << n_rejected_nmvtx
              << ", rejected nINTT=" << n_rejected_nintt
              << ", rejected vertex multiplicity="
              << n_rejected_vertex_multiplicity
              << ", rejected vertex high-pT multiplicity="
              << n_rejected_vertex_high_pt_multiplicity
              << ", missing matching vertices="
              << n_missing_matching_vertices
              << ", tracks with absent vertex="
              << n_track_vertex_absent
              << ", tracks with vertex crossing mismatch="
              << n_track_vertex_crossing_mismatch
              << ", absent-vertex tracks with unassigned vertex ID="
              << n_absent_vertex_id_unassigned
              << ", crossing difference underflow (< -50)="
              << n_crossing_diff_underflow
              << ", crossing difference overflow (> 50)="
              << n_crossing_diff_overflow
              << ", rejected vertex position="
              << n_rejected_vertex_position
              << ", insufficient states=" << n_insufficient_track_states
              << ", propagation failures=" << n_propagation_failures
              << ", non-finite vertex residuals="
              << n_nonfinite_vertex_residuals
              << ", unavailable track_ntp tree="
              << n_unavailable_track_ntuple
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::End(PHCompositeNode* /*unused*/)
{
  delete _mille;
  m_constraintFile.close();

  if (m_file)
  {
    if (m_ntuple)
    {
      m_ntuple->Write();
    }
    if (track_ntp)
    {
      track_ntp->Write();
    }
    m_file->Write();
    m_file->Close();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeMilleFiles::GetNodes(PHCompositeNode* topNode)
{
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name);
  if (!_track_map)
  {
    std::cout
        << PHWHERE << "\n"
        << "\tCan't find track map:\n"
        << "\t" << m_track_map_name << "\n"
        << "\tAborting\n"
        << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (m_useEventVertex && !m_vertex_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node SvtxVertexMap" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tGeometry)
  {
    std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _state_map = findNode::getClass<SvtxAlignmentStateMap>(topNode, m_state_map_name);
  if (!_state_map)
  {
    std::cout
        << PHWHERE << "\n"
        << "\tCan't find state map:\n"
        << "\t" << m_state_map_name << "\n"
        << "\tAborting\n"
        << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool MakeMilleFiles::getLocalVtxDerivativesXY(SvtxTrack* track,
                                              ActsPropagator& propagator,
                                              const Acts::Vector3& vertex,
                                              float lclvtx_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NLOC],
                                              Acts::BoundSquareMatrix* perigeeCov)
{

  //! Get the first track state beyond the vertex, which will be the
  //! innermost track state and propagate it to the vertex surface to
  //! get the jacobian at the vertex
  auto* firststate = (*std::next(track->begin_states(), 1)).second;

  const auto ckey = firststate->get_cluskey();
  auto *const cluster = _cluster_map->findCluster(ckey);
  const auto surf = _tGeometry->maps().getSurface(ckey, cluster);

  const auto param = propagator.makeTrackParams(firststate, track->get_charge(), surf).value();
  const auto perigee = propagator.makeVertexSurface(vertex);
  const auto actspropagator = propagator.makePropagator();
  ActsPropagator::SphenixPropagator::Options options(
      _tGeometry->geometry().getGeoContext(),
      _tGeometry->geometry().magFieldContext);
  const auto intersection = perigee->intersect(
      _tGeometry->geometry().getGeoContext(),
      param.position(_tGeometry->geometry().getGeoContext()),
      param.momentum(), Acts::BoundaryTolerance::None(),
      0.1 * Acts::UnitConstants::mm).closest();
  if (!intersection.isValid())
  {
    if (Verbosity() > 1)
    {
      std::cout << "MakeMilleFiles::getLocalVtxDerivativesXY - no valid "
                << "intersection with the vertex perigee surface for track "
                << track->get_id() << std::endl;
    }
    return false;
  }
  options.direction =
      Acts::Direction::fromScalarZeroAsPositive(intersection.pathLength());

  const auto result = actspropagator.propagate(param, *perigee, options);

  if (result.ok())
  {
    //auto jacobian = *result.value().transportJacobian;

    // Bound covariance (Acts units: mm, rad, GeV) of the track on the vertex
    // perigee surface, used as the track covariance for the pseudo-measurements.
    // Left as NaN if unavailable so the pseudo-measurement step rejects it.
    if (perigeeCov)
    {
      const auto& endParams = result.value().endParameters;
      if (endParams && endParams->covariance())
      {
        *perigeeCov = *endParams->covariance();
      }
      else
      {
        perigeeCov->setConstant(std::numeric_limits<double>::quiet_NaN());
      }
    }

    Eigen::Matrix<double, 2, 6> projector = Eigen::Matrix<double, 2, 6>::Zero();
    projector(0, 0) = 1;
    projector(1, 1) = 1;
    Eigen::Matrix<double, 2, 6> deriv = projector;// * jacobian;
    if (Verbosity() > 2)
    {
      std::cout << "local vtxderiv (full, before truncation) " << std::endl
                << deriv << std::endl;
    }
    // Temporary: the cluster rows only use local labels 1-2 (loc0, loc1), so
    // keep only d(d0, z0)/d(loc0, loc1) here. Columns for phi, theta, q/p and t
    // are not constrained by any cluster row and make pede's local matrix
    // rank deficient.
    deriv.rightCols<4>().setZero();

    for (int i = 0; i < deriv.rows(); i++)
    {
      for (int j = 0; j < deriv.cols(); j++)
      {
        lclvtx_derivative[i][j] = deriv(i, j);
      }
    }

    return true;
  }

  return false;
}

void MakeMilleFiles::getGlobalVtxDerivativesXY(SvtxTrack* track,
                                               const Acts::Vector3& vertex,
                                               float glblvtx_derivative[SvtxAlignmentState::NRES][3])
{
  Acts::SquareMatrix3 identity = Acts::SquareMatrix3::Identity();

  Acts::Vector3 projx = Acts::Vector3::Zero();
  Acts::Vector3 projy = Acts::Vector3::Zero();
  getProjectionVtxXY(track, vertex, projx, projy);

  glblvtx_derivative[0][0] = identity.col(0).dot(projx);
  glblvtx_derivative[0][1] = identity.col(1).dot(projx);
  glblvtx_derivative[0][2] = identity.col(2).dot(projx);
  glblvtx_derivative[1][0] = identity.col(0).dot(projy);
  glblvtx_derivative[1][1] = identity.col(1).dot(projy);
  glblvtx_derivative[1][2] = identity.col(2).dot(projy);

  /// swap sign to fit expectation of pede to have derivative of fit rather than
  /// residual
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      glblvtx_derivative[i][j] *= -1;
    }
  }
}

bool MakeMilleFiles::getLocalVtxDerivativesXYpsuedo(SvtxTrack* track,
                                                    const SvtxVertex* vertex,
                                                    const Acts::BoundSquareMatrix& perigeeCov,
                                                    float localderpsuedo[SvtxAlignmentState::NLOC][SvtxAlignmentState::NLOC],
                                                    float localmeaserrpsuedo[SvtxAlignmentState::NLOC])
{
  using LocalCov = Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>;
  using LocalVec = Eigen::Matrix<double, SvtxAlignmentState::NLOC, 1>;
  //using ResCov = Eigen::Matrix<double, SvtxAlignmentState::NRES, SvtxAlignmentState::NRES>;
  //using ResCov2 = Eigen::Matrix<double, 3,3>;

  for (int i = 0; i < SvtxAlignmentState::NLOC; ++i)
  {
    localmeaserrpsuedo[i] = 0;
    for (int j = 0; j < SvtxAlignmentState::NLOC; ++j)
    {
      localderpsuedo[i][j] = 0;
    }
  }

  if (!track || !vertex)
  {
    std::cout << "MakeMilleFiles::getLocalVtxDerivativesXYpsuedo - missing "
              << (track ? "vertex" : "track") << std::endl;
    return false;
  }
  if (!perigeeCov.allFinite())
  {
    std::cout << "MakeMilleFiles::getLocalVtxDerivativesXYpsuedo - no finite "
              << "perigee covariance for track " << track->get_id() << std::endl;
    return false;
  }

  //! Track covariance of the bound parameters on the vertex perigee surface.
  //! This plays the role of state.covariance() for a cluster state.
  const LocalCov updatedCovariance =
      perigeeCov.topLeftCorner<SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>();

  //! Vertex covariance expressed in the perigee (d0, z0) frame.
  //! A vertex shift dv changes the perigee parameters by -P dv, with
  //!   row 0: Acts loc0 axis, (z x p)/|z x p|
  //!   row 1: z - cot(theta) n, n = unit transverse momentum
  //! Row 1 carries the dip correction: moving the vertex along the transverse
  //! track direction moves the point of closest approach, and hence z0.
  //! These are the projections of getProjectionVtxXY (ATLAS eq. 31), with
  //! row 0 sign-matched to the Acts perigee frame so the d0-z0 correlation
  //! has the same sign convention as the track covariance.
  //const Acts::Vector3 direction =
  //    Acts::Vector3(track->get_px(), track->get_py(), track->get_pz()).normalized();
  //const Acts::Vector3 zhat = Acts::Vector3::UnitZ();
  //const Acts::Vector3 transverse = Acts::Vector3(direction.x(), direction.y(), 0).normalized();
  //const double cotTheta = direction.z() / std::hypot(direction.x(), direction.y());

  //Eigen::Matrix<double, SvtxAlignmentState::NRES, 3> vertexProjection;
  //vertexProjection.row(0) = zhat.cross(direction).normalized().transpose();
  //vertexProjection.row(1) = (zhat - cotTheta * transverse).transpose();

  ////! SvtxVertex errors are in cm^2; the Acts track covariance is in mm^2
  //Acts::SquareMatrix3 vertexCov;
  //for (int i = 0; i < 3; ++i)
  //{
  //  for (int j = 0; j < 3; ++j)
  //  {
  //    vertexCov(i, j) = vertex->get_error(i, j) *
  //                      Acts::UnitConstants::cm * Acts::UnitConstants::cm;
  //  }
  //}

  //const ResCov rotatedVertexCov =
  //    vertexProjection * vertexCov * vertexProjection.transpose();
  //
  //Eigen::Matrix<double, 3, 3> vertexProjection2;
  //Acts::Vector3 r = direction.cross(Acts::Vector3(0., 0., 1.));
  //float phi = atan2(r(1), r(0));
  //phi *= -1;
  //vertexProjection2(0, 0) = std::cos(phi);
  //vertexProjection2(0, 1) = -std::sin(phi);
  //vertexProjection2(0, 2) = 0;
  //vertexProjection2(1, 0) = std::sin(phi);
  //vertexProjection2(1, 1) = std::cos(phi);
  //vertexProjection2(1, 2) = 0;
  //vertexProjection2(2, 0) = 0;
  //vertexProjection2(2, 1) = 0;
  //vertexProjection2(2, 2) = 1;
  //const ResCov2 rotatedVertexCov2 = vertexProjection2 * vertexCov * vertexProjection2.transpose();
   
  //std::cout << "MakeMilleFiles::getLocalVtxDerivativesXYpsuedo - track "
  //           << track->get_id() << ", vertex " << vertex->get_id() << std::endl
  //           << "  global vertex covariance (mm^2, x/y/z):" << std::endl
  //           << vertexCov << std::endl
  //           << "  rotated vertex covariance (mm^2, d0/z0 perigee frame):" << std::endl
  //           << rotatedVertexCov << std::endl
  //           << "  rotated vertex covariance 2 (mm^2, phi rotation frame):" << std::endl
  //           << rotatedVertexCov2 << std::endl;

  
  //if (!rotatedVertexCov.allFinite() || rotatedVertexCov.determinant() <= 0)
  //{
  //  std::cout << "MakeMilleFiles::getLocalVtxDerivativesXYpsuedo - WARNING: "
  //            << "rotated vertex covariance is not positive definite" << std::endl;
  //}

  //! The track fit does not include the vertex, so the perigee covariance
  //! already excludes it: use C^-1 directly rather than the complement
  //! C^-1 - H^T V^-1 H used for cluster states.
  const LocalCov regularisedCov =
      regulariseCovariance(updatedCovariance, -1, -1, 1.e-9);

  LocalCov correlationTerm =
      regularisedCov.ldlt().solve(LocalCov::Identity());
  // symmetrise to correct for floating point effects
  correlationTerm = 0.5 * (correlationTerm + correlationTerm.transpose());

  Eigen::SelfAdjointEigenSolver<LocalCov> eigenSolver(correlationTerm);
  if (eigenSolver.info() != Eigen::Success)
  {
    std::cout << " FAILED to find decompose correlation term" << std::endl;
    return false;
  }
  const LocalVec& eigenVals = eigenSolver.eigenvalues();
  const LocalCov& eigenVecs = eigenSolver.eigenvectors();

  /// convert each EV to a pseudo-measurement. Laid out for Mille: row iMeas
  /// is the local derivative of pseudo-measurement iMeas (its eigenvector),
  /// and localmeaserrpsuedo[iMeas] its sigma. Skipped (non-positive) EVs keep
  /// a sigma of 0 and must not be passed to Mille.
  for (int iMeas = 0; iMeas < SvtxAlignmentState::NLOC; ++iMeas)
  {
    // skip negative EV
    if (eigenVals(iMeas) <= 0)
    {
      continue;
    }
    localmeaserrpsuedo[iMeas] = static_cast<float>(1. / std::sqrt(eigenVals(iMeas)));
    for (int iPar = 0; iPar < SvtxAlignmentState::NLOC; ++iPar)
    {
      localderpsuedo[iMeas][iPar] = static_cast<float>(eigenVecs(iPar, iMeas));
    }
  }

  if (Verbosity() > 2)
  {
    std::cout << "  vertex pseudo-measurement eigenvalues: "
              << eigenVals.transpose() << std::endl
              << "  vertex pseudo-measurement (sigma: derivatives):" << std::endl;
    for (int iMeas = 0; iMeas < SvtxAlignmentState::NLOC; ++iMeas)
    {
      std::cout << "    " << localmeaserrpsuedo[iMeas] << ":";
      for (int iPar = 0; iPar < SvtxAlignmentState::NLOC; ++iPar)
      {
        std::cout << " " << localderpsuedo[iMeas][iPar];
      }
      std::cout << std::endl;
    }
  }

  return true;
}

std::pair<float,float> MakeMilleFiles::getRotatedVtxCovarianceDiagonalPerigee(SvtxTrack* track,
                                                    const SvtxVertex* vertex)
{
  using ResCov = Eigen::Matrix<double, SvtxAlignmentState::NRES, SvtxAlignmentState::NRES>;


  //! Vertex covariance expressed in the perigee (d0, z0) frame.
  //! A vertex shift dv changes the perigee parameters by -P dv, with
  //!   row 0: Acts loc0 axis, (z x p)/|z x p|
  //!   row 1: z - cot(theta) n, n = unit transverse momentum
  //! Row 1 carries the dip correction: moving the vertex along the transverse
  //! track direction moves the point of closest approach, and hence z0.
  //! These are the projections of getProjectionVtxXY (ATLAS eq. 31), with
  //! row 0 sign-matched to the Acts perigee frame so the d0-z0 correlation
  //! has the same sign convention as the track covariance.
  const Acts::Vector3 direction =
      Acts::Vector3(track->get_px(), track->get_py(), track->get_pz()).normalized();
  const Acts::Vector3 zhat = Acts::Vector3::UnitZ();
  const Acts::Vector3 transverse = Acts::Vector3(direction.x(), direction.y(), 0).normalized();
  const double cotTheta = direction.z() / std::hypot(direction.x(), direction.y());

  Eigen::Matrix<double, SvtxAlignmentState::NRES, 3> vertexProjection;
  vertexProjection.row(0) = zhat.cross(direction).normalized().transpose();
  vertexProjection.row(1) = (zhat - cotTheta * transverse).transpose();

  //! SvtxVertex errors are in cm^2; the Acts track covariance is in mm^2
  Acts::SquareMatrix3 vertexCov;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      vertexCov(i, j) = vertex->get_error(i, j) *
                        Acts::UnitConstants::cm * Acts::UnitConstants::cm;
    }
  }

  ResCov rotatedVertexCov =
      vertexProjection * vertexCov * vertexProjection.transpose();

  return std::make_pair(pow(rotatedVertexCov(0, 0),.5), pow(rotatedVertexCov(1, 1),.5));
}
Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> MakeMilleFiles::regulariseCovariance(const Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>& inputCov,
                                                 double conditionCutOff,
                                                 double removeHugeLeading,
                                                 double stabilisationDiag)
{
  Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> out =
      Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>::Zero(inputCov.rows(), inputCov.cols());

  /// optionally add a tiny diagonal matrix for additional stabilisation
  Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> regularisation =
      Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>::Zero(inputCov.rows(), inputCov.cols());
  if (stabilisationDiag > 0) {
    regularisation = stabilisationDiag * Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>::Identity(
                                             inputCov.rows(), inputCov.cols());
  }
  auto eigensolver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>>(
      inputCov + regularisation);

  if (eigensolver.info() != Eigen::Success) {
    std::cout << " FAILED to find eigenvec" << std::endl;
    return out;
  }
  auto eigenVals = eigensolver.eigenvalues();
  auto eigenVecs = eigensolver.eigenvectors();

  std::size_t maxIndex = eigenVals.size() - 1;
  double lambdaMax = eigenVals(eigenVals.size() - 1);
  // check for a huge leading eigenvalue - this happens when the time coordinate
  // is unconstrained
  if (removeHugeLeading > 0 &&
      lambdaMax > removeHugeLeading * eigenVals(maxIndex - 1)) {
    --maxIndex;
    lambdaMax = eigenVals(maxIndex);
  }
  double lambdaMin = eigenVals(0);
  if (conditionCutOff > 0) {
    lambdaMin = conditionCutOff * lambdaMax;
  }
  // clamp the EV to the permitted interval
  for (Eigen::Index i = 0; i < eigenVals.size(); ++i) {
    out(i, i) = std::clamp(eigenVals(i), lambdaMin, lambdaMax);
  }
  // then return the regularised matrix as V D V^T
  out = out * eigenVecs.transpose();
  out = eigenVecs * out;
  return out;
}

Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> MakeMilleFiles::getInverseComplement(Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>& target,
                                           const Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>& existing_sol)
{
  Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> Rhs =
      Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC>::Identity(target.rows(), target.cols()) -
      target * existing_sol;
  // call decomposition from Eigen (prefer over llt for semi-def. matrices)
  auto LDLT = target.ldlt();
  Eigen::Matrix<double, SvtxAlignmentState::NLOC, SvtxAlignmentState::NLOC> solution = LDLT.solve(Rhs);
  // finally, symmetrise to correct for floating point effects
  solution = 0.5 * (solution + solution.transpose());
  return solution;
}

void MakeMilleFiles::getProjectionVtxXY(SvtxTrack* track,
                                        const Acts::Vector3& vertex,
                                        Acts::Vector3& projx,
                                        Acts::Vector3& projy)
{
  Acts::Vector3 tangent(track->get_px(), track->get_py(), track->get_pz());
  Acts::Vector3 normal(track->get_px(), track->get_py(), 0);

  tangent /= tangent.norm();
  normal /= normal.norm();

  Acts::Vector3 localx(1, 0, 0);
  Acts::Vector3 localz(0, 0, 1);

  Acts::Vector3 xglob = localToGlobalVertex(track, vertex, localx);
  Acts::Vector3 yglob = localz + vertex;

  Acts::Vector3 X = (xglob - vertex) / (xglob - vertex).norm();
  Acts::Vector3 Y = (yglob - vertex) / (yglob - vertex).norm();

  // see equation 31 of the ATLAS paper (and discussion) for this
  projx = X - (tangent.dot(X) / tangent.dot(normal)) * normal;
  projy = Y - (tangent.dot(Y) / tangent.dot(normal)) * normal;

  return;
}
Acts::Vector3 MakeMilleFiles::localToGlobalVertex(SvtxTrack* track,
                                                  const Acts::Vector3& vertex,
                                                  const Acts::Vector3& localx) const
{
  Acts::Vector3 mom(track->get_px(),
                    track->get_py(),
                    track->get_pz());

  Acts::Vector3 r = mom.cross(Acts::Vector3(0., 0., 1.));
  float phi = atan2(r(1), r(0));
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;

  rot(0, 0) = std::cos(phi);
  rot(0, 1) = -std::sin(phi);
  rot(0, 2) = 0;
  rot(1, 0) = std::sin(phi);
  rot(1, 1) = std::cos(phi);
  rot(1, 2) = 0;
  rot(2, 0) = 0;
  rot(2, 1) = 0;
  rot(2, 2) = 1;

  rot_T = rot.transpose();

  Acts::Vector3 pos_R = rot * localx;

  pos_R += vertex;

  return pos_R;
}
Acts::Vector3 MakeMilleFiles::getEventVertex(const SvtxTrack* track, const SvtxVertex*& svtxVertexOut) const
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  svtxVertexOut = nullptr;

  if (!track || !m_vertex_map)
  {
    return Acts::Vector3::Constant(nan);
  }

  const SvtxVertex* vertex =
      m_vertex_map->get(track->get_vertex_id());
  if (!vertex ||
      vertex->get_beam_crossing() != track->get_crossing())
  {
    return Acts::Vector3::Constant(nan);
  }
  svtxVertexOut = vertex;
  return Acts::Vector3(
      svtxVertexOut->get_x(),
      svtxVertexOut->get_y(),
      svtxVertexOut->get_z());
}

void MakeMilleFiles::addTrackToMilleFile(SvtxAlignmentStateMap::StateVec& statevec)
{
  for (auto *state : statevec)
  {
    TrkrDefs::cluskey ckey = state->get_cluster_key();

    if (Verbosity() > 2)
    {
      std::cout << "adding state for ckey " << ckey << " with hitsetkey "
                << (int) TrkrDefs::getHitSetKeyFromClusKey(ckey) << std::endl;
    }
    // The global alignment parameters are given initial values of zero by default, we do not specify them
    // We identify the global alignment parameters for this surface

    TrkrCluster* cluster = _cluster_map->findCluster(ckey);
    const unsigned int layer = TrkrDefs::getLayer(ckey);
    const unsigned int trkrid = TrkrDefs::getTrkrId(ckey);
     uint8_t stave = -1;
      uint8_t chip = -1;
      if (trkrid == TrkrDefs::mvtxId)
      {
        // need stave to get clamshell
         stave = MvtxDefs::getStaveId(ckey);
         chip = MvtxDefs::getChipId(ckey);
      }
      else if(trkrid == TrkrDefs::inttId)
      {
        stave = InttDefs::getLadderPhiId(ckey);
        chip = InttDefs::getLadderZId(ckey);
      }
    if(m_ignore_tpc && trkrid == TrkrDefs::tpcId) {
      continue;
}
    const SvtxAlignmentState::ResidualVector residual = state->get_residual();// / Acts::UnitConstants::cm;
    //auto acts_pars = state->parameters();
    //std::cout<<"acts_par(0): "<<acts_pars(0)<<std::endl;
    const Acts::Vector3 global = _tGeometry->getGlobalPosition(ckey, cluster);

    // need standard deviation of measurements
    SvtxAlignmentState::ResidualVector clus_sigma = SvtxAlignmentState::ResidualVector::Zero();

    //double clusRadius = sqrt(global[0] * global[0] + global[1] * global[1]);
    double clusphi = atan2(global[1] , global[0]);
    //auto para_errors = _ClusErrPara.get_clusterv5_modified_error(cluster, clusRadius, ckey);
    //double phierror = sqrt(para_errors.first);
    //double zerror = sqrt(para_errors.second);
    double phierror = cluster->getRPhiError();
    double zerror = cluster->getZError();
    clus_sigma(1) = zerror * Acts::UnitConstants::cm;
    clus_sigma(0) = phierror * Acts::UnitConstants::cm;

    if (std::isnan(clus_sigma(0)) ||
        std::isnan(clus_sigma(1)))
    {
      continue;
    }

    auto surf = _tGeometry->maps().getSurface(ckey, cluster);

    int glbl_label[SvtxAlignmentState::NGL];
    if (layer < 3)
    {
      AlignmentDefs::getMvtxGlobalLabels(surf, ckey, glbl_label, mvtx_group);
    }
    else if (layer > 2 && layer < 7)
    {
      AlignmentDefs::getInttGlobalLabels(surf, ckey, glbl_label, intt_group);
    }
    else if (layer < 55)
    {
      AlignmentDefs::getTpcGlobalLabels(surf, ckey, glbl_label, tpc_group);
    }
    else if (layer < 57)
    {
      AlignmentDefs::getMMGlobalLabels(surf, glbl_label, mms_group);
    }

    if (Verbosity() > 1)
    {
      std::cout << std::endl;
    }

    float glbl_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NGL]{};
    float lcl_derivative[SvtxAlignmentState::NRES][SvtxAlignmentState::NLOC]{};
    float lcl_derivative_psuedo[SvtxAlignmentState::NLOC][SvtxAlignmentState::NLOC]{};
    float lcl_meas_psuedo[SvtxAlignmentState::NLOC]{};
    float glbl_derivative_dummy_empty[SvtxAlignmentState::NGL]{};
    int glbl_label_dummy_empty[SvtxAlignmentState::NGL]{};
    
    /// For N residual local coordinates x, z
    for (int i = 0; i < SvtxAlignmentState::NRES; ++i)
    {
     
      // Add the measurement separately for each coordinate direction to Mille
      for (int j = 0; j < SvtxAlignmentState::NGL; ++j)
      {
        glbl_derivative[i][j] = state->get_global_derivative_matrix()(i, j);
        if (is_layer_fixed(layer) ||
            is_layer_param_fixed(layer, j, fixed_layer_gparams))
        {
          glbl_derivative[i][j] = 0.;
        }

        if (trkrid == TrkrDefs::mvtxId)
        {
          // need stave to get clamshell
           auto stave4clam = MvtxDefs::getStaveId(ckey);
          auto clamshell = AlignmentDefs::getMvtxClamshell(layer, stave4clam);
          if (is_layer_param_fixed(layer, j, fixed_layer_gparams) ||
              is_mvtx_layer_fixed(layer, clamshell))
          {
            glbl_derivative[i][j] = 0;
          }
        }
        else if (trkrid == TrkrDefs::inttId)
        {
          
          if (is_layer_param_fixed(layer, j, fixed_layer_gparams) ||
              is_layer_fixed(layer))
          {
            glbl_derivative[i][j] = 0;
          }
        }
        if (trkrid == TrkrDefs::tpcId)
        {
          auto sector = TpcDefs::getSectorId(ckey);
          auto side = TpcDefs::getSide(ckey);
          if (is_layer_param_fixed(layer, j, fixed_layer_gparams) ||
              is_tpc_sector_fixed(layer, sector, side) ||
              is_layer_fixed(layer))
          {
            glbl_derivative[i][j] = 0.0;
          }
        }
      }

      for (int j = 0; j < SvtxAlignmentState::NLOC; ++j)
      {
        lcl_derivative[i][j] = state->get_local_derivative_matrix()(i, j);

        if (is_layer_param_fixed(layer, j, fixed_layer_lparams))
        {
          lcl_derivative[i][j] = 0.;
        }
      }
      if (Verbosity() > 2)
      {
        std::cout << "coordinate " << i << " has residual " << residual(i) << " and clus_sigma " << clus_sigma(i) << std::endl
                  << "global deriv " << std::endl;

        for (float k : glbl_derivative[i])
        {
          if (k > 0 || k < 0)
          {
            std::cout << "NONZERO GLOBAL DERIVATIVE" << std::endl;
          }
          std::cout << k << ", ";
        }
        std::cout << std::endl
                  << "local deriv " << std::endl;
        for (float k : lcl_derivative[i])
        {
          std::cout << k << ", ";
        }
        std::cout << std::endl;
      }

      if (clus_sigma(i) < 1.0 &&  // discards crazy clusters
          std::isfinite(residual(i)) &&
          allFinite(lcl_derivative[i], SvtxAlignmentState::NLOC) &&
          allFinite(glbl_derivative[i], SvtxAlignmentState::NGL))
      {
        if (Verbosity() > 3)
        {
          std::cout << "ckey " << ckey << " and layer " << layer << " buffers:" << std::endl;
          AlignmentDefs::printBuffers(i, residual, clus_sigma, lcl_derivative[i], glbl_derivative[i], glbl_label);
        }
        float errinf = 1.0;
        if (m_layerMisalignment.contains(layer))
        {
          errinf = m_layerMisalignment.find(layer)->second;
        }

        _mille->mille(SvtxAlignmentState::NLOC, lcl_derivative[i], SvtxAlignmentState::NGL, glbl_derivative[i], glbl_label, residual(i), errinf * clus_sigma(i));
      }
    }

    // For 6 local track parameters, the associated psuedo mille measurements that complete the local derivative
    for (int i = 0; i < SvtxAlignmentState::NLOC; ++i)
    {
      for (int j = 0; j < SvtxAlignmentState::NLOC; ++j)
      {
        lcl_derivative_psuedo[i][j] = state->get_local_derivative_psuedo_matrix()(i, j);
      }
      lcl_meas_psuedo[i] = state->get_local_psuedo_measurement_err()(i);

      if (std::isfinite(lcl_meas_psuedo[i]) &&
          allFinite(lcl_derivative_psuedo[i], SvtxAlignmentState::NLOC))
      {
        _mille->mille(SvtxAlignmentState::NLOC, lcl_derivative_psuedo[i], SvtxAlignmentState::NGL, glbl_derivative_dummy_empty, glbl_label_dummy_empty, 0.0, lcl_meas_psuedo[i]);
      }
    }

//    // DEBUG-PSEUDO-TEST begin: pseudo-measurement sigmas as read by MakeMilleFiles
//    {
//      static int n_debug_pseudo_reads = 0;
//      if (n_debug_pseudo_reads < 20)
//      {
//        ++n_debug_pseudo_reads;
//        std::cout << "DEBUG-PSEUDO MakeMilleFiles ckey " << ckey << " read sigma:";
//        for (const float sigma : lcl_meas_psuedo)
//        {
//          std::cout << " " << sigma;
//        }
//        std::cout << " (Mille::mille drops rows with sigma <= 0)" << std::endl;
//      }
//    }
//    // DEBUG-PSEUDO-TEST end

    float ntp_data[] = {
      (float) layer, (float) stave, (float) chip, (float) residual(0), (float) residual(1), (float) clusphi, (float) global[0], (float) global[1], (float) global[2],
      (float) clus_sigma(0), (float) clus_sigma(1),
      lcl_derivative[0][0], lcl_derivative[0][1], lcl_derivative[0][2], lcl_derivative[0][3], lcl_derivative[0][4], lcl_derivative[0][5],
      glbl_derivative[0][0], glbl_derivative[0][1], glbl_derivative[0][2], glbl_derivative[0][3], glbl_derivative[0][4], glbl_derivative[0][5],
      lcl_derivative[1][0], lcl_derivative[1][1], lcl_derivative[1][2], lcl_derivative[1][3], lcl_derivative[1][4], lcl_derivative[1][5],
      glbl_derivative[1][0], glbl_derivative[1][1], glbl_derivative[1][2], glbl_derivative[1][3], glbl_derivative[1][4], glbl_derivative[1][5],
    };

    if (m_ntuple) { m_ntuple->Fill(ntp_data);
}
  }

  return;
}

bool MakeMilleFiles::is_layer_fixed(unsigned int layer)
{
  bool ret = false;
  auto it = fixed_layers.find(layer);
  if (it != fixed_layers.end())
  {
    ret = true;
  }

  return ret;
}
void MakeMilleFiles::set_layers_fixed(unsigned int minlayer,
                                      unsigned int maxlayer)
{
  for (unsigned int i = minlayer; i < maxlayer; i++)
  {
    fixed_layers.insert(i);
  }
}
void MakeMilleFiles::set_layer_fixed(unsigned int layer)
{
  fixed_layers.insert(layer);
}

bool MakeMilleFiles::is_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell)
{
  bool ret = false;

  std::pair pair = std::make_pair(layer, clamshell);
  auto it = fixed_mvtx_layers.find(pair);
  if (it != fixed_mvtx_layers.end())
  {
    ret = true;
  }

  return ret;
}
void MakeMilleFiles::set_mvtx_layer_fixed(unsigned int layer, unsigned int clamshell)
{
  fixed_mvtx_layers.insert(std::make_pair(layer, clamshell));
}
bool MakeMilleFiles::is_layer_param_fixed(unsigned int layer, unsigned int param, std::set<std::pair<unsigned int, unsigned int>>& param_fixed)
{
  bool ret = false;
  std::pair<unsigned int, unsigned int> pair = std::make_pair(layer, param);
  auto it = param_fixed.find(pair);
  if (it != param_fixed.end())
  {
    ret = true;
  }

  return ret;
}

void MakeMilleFiles::set_layer_gparam_fixed(unsigned int layer, unsigned int param)
{
  fixed_layer_gparams.insert(std::make_pair(layer, param));
}
void MakeMilleFiles::set_layer_lparam_fixed(unsigned int layer, unsigned int param)
{
  fixed_layer_lparams.insert(std::make_pair(layer, param));
}
bool MakeMilleFiles::is_tpc_sector_fixed(unsigned int layer, unsigned int sector, unsigned int side)
{
  bool ret = false;
  unsigned int region = AlignmentDefs::getTpcRegion(layer);
  unsigned int subsector = region * 24 + side * 12 + sector;
  auto it = fixed_sectors.find(subsector);
  if (it != fixed_sectors.end())
  {
    ret = true;
  }

  return ret;
}
