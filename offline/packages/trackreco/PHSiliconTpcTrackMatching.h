// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONTPCTRACKMATCHING_H
#define PHSILICONTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>
#include <trackbase/ActsGeometry.h>

#include <map>
#include <string>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class TrkrClusterContainer;
class TF1;
class TrkrClusterCrossingAssoc;

class PHSiliconTpcTrackMatching : public SubsysReco, public PHParameterInterface
{
 public:
  PHSiliconTpcTrackMatching(const std::string &name = "PHSiliconTpcTrackMatching");

  ~PHSiliconTpcTrackMatching() override;

  void SetDefaultParameters() override;

  void set_phi_search_window(const double win) { _phi_search_win = win; }
  void set_eta_search_window(const double win) { _eta_search_win = win; }
  void set_x_search_window(const double win) { _x_search_win = win; }
  void set_y_search_window(const double win) { _y_search_win = win; }
  void set_z_search_window(const double win) { _z_search_win = win; }

  void set_match_window_function_pars(const double a, const double b, const double ptmin)
  {
    _match_function_a = a;
    _match_function_b = b;
    _match_function_ptmin = ptmin;
  }
  void set_use_old_matching(const bool flag) { _use_old_matching = flag; }

  void set_test_windows_printout(const bool test) { _test_windows = test; }
  void set_pp_mode(const bool flag) { _pp_mode = flag; }
  void set_use_intt_crossing(const bool flag) { _use_intt_crossing = flag; }

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  void fieldMap(std::string &fieldmap) { m_fieldMap = fieldmap; }

  void set_silicon_track_map_name(const std::string &map_name) { _silicon_track_map_name = map_name; }
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void SetIteration(int iter) { _n_iteration = iter; }

 private:
  int GetNodes(PHCompositeNode *topNode);

  void findEtaPhiMatches(std::set<unsigned int> &tpc_matched_set,
                         std::set<unsigned int> &tpc_unmatched_set,
                         std::multimap<unsigned int, unsigned int> &tpc_matches);
  std::vector<short int> getInttCrossings(TrackSeed *si_track);
  void checkCrossingMatches(std::multimap<unsigned int, unsigned int> &tpc_matches);
  short int getCrossingIntt(TrackSeed *_tracklet_si);
  // void findCrossingGeometrically(std::multimap<unsigned int, unsigned int> tpc_matches);
  short int findCrossingGeometrically(unsigned int tpc_id, unsigned int si_id);
  double getBunchCrossing(unsigned int trid, double z_mismatch);
  double getMatchingInflationFactor(double tpc_pt);

  // default values, can be replaced from the macro
  double _phi_search_win = 0.01;
  double _eta_search_win = 0.004;
  double _x_search_win = 0.3;
  double _y_search_win = 0.3;
  double _z_search_win = 0.4;

  double _match_function_a = 1.0;
  double _match_function_b = 5.0;
  double _match_function_pow = 1.0;
  double _match_function_ptmin = 0.15;
  bool _use_old_matching = false;  // normally false

  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_track_map{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrackSeed *_tracklet_tpc{nullptr};
  TrackSeed *_tracklet_si{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};
  TrkrClusterCrossingAssoc *_cluster_crossing_map{nullptr};

  std::map<unsigned int, double> _z_mismatch_map;

  //  double _collision_rate = 50e3;  // input rate for phi correction
  //  double _reference_collision_rate = 50e3;  // reference rate for phi correction
  //  double _si_vertex_dzmax = 0.25;  // mm
  double crossing_period = 106.0;  // ns
  double fieldstrength{std::numeric_limits<double>::quiet_NaN()};

  bool _test_windows = false;
  bool _pp_mode = false;
  bool _use_intt_crossing = true;  // should always be true except for testing

  int _n_iteration = 0;
  std::string _track_map_name = "TpcTrackSeedContainer";
  std::string _silicon_track_map_name = "SiliconTrackSeedContainer";
  std::string m_fieldMap = "1.4";
};

#endif  //  PHSILICONTPCTRACKMATCHING_H
