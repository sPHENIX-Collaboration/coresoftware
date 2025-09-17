// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKFITTINGQA_H
#define TRACKFITTINGQA_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <limits>

class TH1;
class PHCompositeNode;

class TrackFittingQA : public SubsysReco
{
 public:
  TrackFittingQA(const std::string& name = "TrackFittingQA");
  ~TrackFittingQA() override = default;

  /// Histogram ranges (see member declarations for default values)
  void set_quality_range ( float const& quality_range ) { m_quality_range = quality_range; }
  void set_momentum_range ( float const& momentum_range ) { m_momentum_range = momentum_range; }
  //...

  /// Cuts (bounds are inclusive, e.g. set_min_intt_states(2) means 1 state tracks are ignored, while 2 state tracks are kept)
  /// Defaults are all-inclusive
  void set_min_quality ( float const& min_quality ) { m_min_quality = min_quality; }
  void set_max_quality ( float const& max_quality ) { m_max_quality = max_quality; }
  void set_min_p ( float const& min_p ) { m_min_p = min_p; }
  void set_min_pt ( float const& min_pt ) { m_min_pt = min_pt; }
  void set_max_abs_eta ( float const& max_abs_eta ) { m_max_abs_eta = max_abs_eta; }
  void set_min_intt_states ( int const& min_intt_states ) { m_min_intt_states = min_intt_states; }
  void set_min_mvtx_states ( int const& min_mvtx_states ) { m_min_mvtx_states = min_mvtx_states; }
  void set_min_tpc_states ( int const& min_tpc_states ) { m_min_tpc_states = min_tpc_states; }
  void set_min_tpot_states ( int const& min_tpot_states ) { m_min_tpot_states = min_tpot_states; }
  void set_min_crossing ( short const& min_crossing ) { m_min_crossing = min_crossing; }
  void set_max_crossing ( short const& max_crossing ) { m_max_crossing = max_crossing; }

  /// sets the name of node to retrieve the track map from (default member value is "SvtxTrackMap")
  void set_track_map_name(std::string const& track_map_node_name) { m_track_map_node_name = track_map_node_name; }

  /// sets the name of node to retrieve the state map from (default member value is "SvtxAlignmentStateMap")
  //	void set_state_map_name (std::string const& state_map_node_name) {m_state_map_node_name = state_map_node_name;}

  int Init(PHCompositeNode*) override;
  int InitRun(PHCompositeNode*) override;

  int process_event(PHCompositeNode*) override;

  int ResetEvent(PHCompositeNode*) override;

  int EndRun(const int runnumber) override;
  int End(PHCompositeNode*) override;

  int Reset(PHCompositeNode*) override;
  void Print(std::string const& = "ALL") const override;

 private:
  float m_momentum_range = 12.0;
  float m_quality_range = 100.0;

  TH1* m_quality_hist[2]{};
  TH1* m_p_hist[2]{};
  TH1* m_pt_hist[2]{};
  TH1* m_pt_err_hist[2]{};
  TH1* m_eta_hist[2]{};
  TH1* m_phi_eta_hist[2]{};
  TH1* m_intt_states_hist[2]{};
  TH1* m_mvtx_states_hist[2]{};
  TH1* m_tpc_states_hist[2]{};
  TH1* m_tpot_states_hist[2]{};
  // ...
  // residual plots as a function of tpc sector (Mariia)

  std::string m_track_map_node_name = "SvtxTrackMap";

  /// Cuts
  float m_min_quality{0};
  float m_max_quality{std::numeric_limits<float>::max()};
  float m_min_p{0};
  float m_min_pt{0};
  float m_max_abs_eta{std::numeric_limits<float>::max()};
  int m_min_intt_states{0};
  int m_min_mvtx_states{0};
  int m_min_tpc_states{0};
  int m_min_tpot_states{0};
  short m_min_crossing{0};
  short m_max_crossing{std::numeric_limits<short>::max()};
};

#endif  // TRACKFITTINGQA_H
