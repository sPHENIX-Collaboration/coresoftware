// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKRECO_PHMICROMEGASTPCTRACKMATCHING_H
#define TRACKRECO_PHMICROMEGASTPCTRACKMATCHING_H

#include <tpc/TpcGlobalPositionWrapper.h>

#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class ActsGeometry;
class TrkrClusterContainer;
class TrkrClusterIterationMapv1;
class TrackSeedContainer;
class PHCompositeNode;
class PHG4CylinderGeomContainer;
class TrackSeed;
class TrkrCluster;
class TF1;
class TH1;

class PHMicromegasTpcTrackMatching : public SubsysReco
{
 public:
  PHMicromegasTpcTrackMatching(const std::string& name = "PHMicromegasTpcTrackMatching");
  ~PHMicromegasTpcTrackMatching() override = default;

  void set_rphi_search_window_lyr1(const double win) { _rphi_search_win[0] = win; }
  void set_z_search_window_lyr1(const double win) { _z_search_win[0] = win; }
  void set_rphi_search_window_lyr2(const double win) { _rphi_search_win[1] = win; }
  void set_z_search_window_lyr2(const double win) { _z_search_win[1] = win; }
  void set_min_tpc_layer(const unsigned int layer) { _min_tpc_layer = layer; }
  void set_test_windows_printout(const bool test) { _test_windows = test; }
  void set_pp_mode(const bool mode) { _pp_mode = mode; }
  void SetIteration(int iter) { _n_iteration = iter; }

  void zeroField(const bool flag) { _zero_field = flag; }

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  // deprecated calls
  inline void set_sc_calib_mode(const bool) {}
  inline void set_collision_rate(const double) {}

 private:
  //! load nodes relevant for the analysis
  int GetNodes(PHCompositeNode* topNode);

  void copyMicromegasClustersToCorrectedMap();

  //! number of layers in the micromegas
  static constexpr unsigned int _n_mm_layers{2};
  
  bool _use_truth_clusters = false;
  bool _zero_field = false;
  TrkrClusterContainer* _cluster_map{nullptr};
  TrkrClusterContainer* _corrected_cluster_map{nullptr};

  TrackSeedContainer* _svtx_seed_map{nullptr};
  TrackSeedContainer* _tpc_track_map{nullptr};
  TrackSeedContainer* _si_track_map{nullptr};

  //! default rphi search window for each layer
  std::array<double, _n_mm_layers> _rphi_search_win{0.25, 13.0};

  //! default z search window for each layer
  std::array<double, _n_mm_layers> _z_search_win{26.0, 0.25};

  // get the cluster list for zeroField
  std::vector<TrkrDefs::cluskey> getTrackletClusterList(TrackSeed* tracklet);
  // range of TPC layers to use in projection to micromegas
  unsigned int _min_tpc_layer{38};

  /// first micromegas layer
  /** it is reset in ::Setup using actual micromegas geometry */
  unsigned int _min_mm_layer{55};

  //! internal event number
  int _event{-1};

  //! micomegas geometry
  PHG4CylinderGeomContainer* _geomContainerMicromegas{nullptr};
  TrkrClusterIterationMapv1* _iteration_map{nullptr};
  int _n_iteration{0};

  //! acts geometry
  ActsGeometry* _tGeometry{nullptr};

  //! tpc global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  //! true to printout actual residuals for testing
  bool _test_windows{false};

  bool _pp_mode{false};
};

#endif  // TRACKRECO_PHMICROMEGASTPCTRACKMATCHING_H
