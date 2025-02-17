// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONTPCTRACKMATCHING_H
#define PHSILICONTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <trackbase/ActsGeometry.h>

#include <map>
#include <string>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class TrkrClusterContainer;
class TF1;
class TrkrClusterCrossingAssoc;
class TFile;
class TNtuple;

class PHSiliconTpcTrackMatching : public SubsysReco, public PHParameterInterface
{
 public:
  PHSiliconTpcTrackMatching(const std::string &name = "PHSiliconTpcTrackMatching");

  ~PHSiliconTpcTrackMatching() override;

  void SetDefaultParameters() override;

  // legacy parameters
  // The legacy code matched tracks as:
  //     |dX| < window * (a+b/pow(pT,c))
  //   for pT < min_pT
  //   otherwise
  //     |dX| < window
  bool _use_legacy_windowing = true;
  void set_use_legacy_windowing (bool set_par=true) { _use_legacy_windowing=set_par; }

  void set_phi_search_window(const double win) { _phi_search_win = win; }
  void set_eta_search_window(const double win) { _eta_search_win = win; }
  void set_x_search_window(const double win) { _x_search_win = win; }
  void set_y_search_window(const double win) { _y_search_win = win; }
  void set_z_search_window(const double win) { _z_search_win = win; }
  void set_crossing_deltaz_max(const double dz) {_crossing_deltaz_max = dz ;}

  float get_phi_search_window() const { return _phi_search_win; }
  float get_eta_search_window() const { return _eta_search_win; }
  float get_x_search_window() const { return _x_search_win; }
  float get_y_search_window() const { return _y_search_win; }
  float get_z_search_window() const { return _z_search_win; }

  // 2024/01/22 update
  struct WindowMatcher {
    // --- option to use the legacy method ---
    bool use_legacy = false;
    double leg_search_win = 1.; // use if use_legacy == true; set in InitRun
    /* PHSiliconTpcTrackMatching* parent_ptr {nullptr}; */
    void set_use_legacy(double _leg_search_win)
    {
      use_legacy=true;
      leg_search_win=_leg_search_win;
    }

    // --- new method, comparing to a+b*exp(c/pT)
    // Each Arr3D object contains a,b,c in order
    // Up to four curved are needed:
    // - for for positive and negative tracks
    // - only one curve if |dX|<fn_max, or two if like  fn_min<dX<fnmax
    using Arr3D = std::array<double,3>;
    std::string print_fn(const Arr3D&);
    Arr3D posLo { 100, 0, 0 }; // (above a2,b2,c2), 100 for |dX|
    Arr3D posHi { 100, 0, 0 }; // (above a3,b3,c3), 100 for use _use_legacy_windowing
    Arr3D negLo { 100, 0, 0 }; // (above a0,b0,c0), 100 for |dX|
    Arr3D negHi { 100, 0, 0 }; // (above a1,b1,c1), 100 for treat all tracks pos Q

    // efficiency flags set during PHSiliconTpcTrackMatching::InitRun()
    bool fabs_max_posQ  = true;
    bool fabs_max_negQ  = true;
    bool negLo_b0 = true;
    bool negHi_b0 = true;
    bool posLo_b0 = true;
    bool posHi_b0 = true;
    double min_pt_posQ  = 0.25; // only grow function windows down to 150 MeV
    double min_pt_negQ  = 0.25; // only grow function windows down to 150 MeV

    WindowMatcher(
        const Arr3D& _posLo={100,0,0},
        const Arr3D& _posHi={100,0,0},
        const Arr3D& _negLo={100,0,0},
        const Arr3D& _negHi={100,0,0},
        const double _min_pt_posQ=0.25,
        const double _min_pt_negQ=0.25)
      : posLo{_posLo}, posHi{_posHi}, negLo{_negLo}, negHi{_negHi},
      min_pt_posQ{_min_pt_posQ}, min_pt_negQ{_min_pt_negQ} {};

    inline double fn_exp(const Arr3D& arr, const bool& b_is_0, double pT) {
      return (b_is_0 ? arr[0] : arr[0]+arr[1]*exp(arr[2]/pT));
    }

    void init_bools(const std::string& which_window="", const bool print=false);

    bool in_window(bool posQ, const double tpc_pt, const double tpc_X, const double si_X);

    // initialize to fn_lo < deltaX < fn_hi for +Q, and fn_lo < deltaX < fn_hi for -Q

    void reset_fns() {
      posLo={100,0,0};
      posHi={100,0,0};
      negLo={100,0,0};
      negHi={100,0,0};
    };

    // same max for |deltaX| for pos and neg Q
    void set_QoverpT_maxabs    (const Arr3D& _posHi, const double _min_pt=0.25)
    { reset_fns(); posHi=_posHi; min_pt_posQ = _min_pt; };

    // same range for deltaX for pos and neg Q
    void set_QoverpT_range     (const Arr3D& _posLo, const Arr3D& _posHi, const double _min_pt=0.25)
    { reset_fns(); posLo=_posLo; posHi=_posHi; min_pt_posQ=_min_pt; };

    // max for |deltaX| for pos Q
    void set_posQoverpT_maxabs (const Arr3D& _posHi, const double _min_pt=0.25)
    { posLo={100.,0.,0.}; posHi=_posHi; min_pt_posQ = _min_pt; };

    // max for |deltaX| for neg Q
    void set_negQoverpT_maxabs (const Arr3D& _negHi, const double _min_pt=0.25)
    { posLo={100.,0.,0.}; negHi=_negHi; min_pt_negQ = _min_pt; };

    // range for deltaX for pos Q
    void set_posQoverpT_range  (const Arr3D& _posLo, const Arr3D& _posHi, const double _min_pt=0.25)
    { posLo=_posLo; posHi=_posHi; min_pt_posQ = _min_pt; };

    // range for deltaX for neg Q
    void set_negQoverpT_range  (const Arr3D& _negLo, const Arr3D& _negHi, const double _min_pt=0.25)
    { negLo=_negLo; negHi=_negHi; min_pt_negQ = _min_pt; };

  };

  // initialize the window matchers with default values
  WindowMatcher window_dx   { {100.,0.,0.}, {5.3, 0., 0.} };
  WindowMatcher window_dy   { {100.,0.,0.}, {5.2, 0., 0.} };
  WindowMatcher window_dz   { {100.,0.,0.}, {0., 2.6, 0.38}, {100,0.,0.,}, {0., 1.45, 0.49} };
  WindowMatcher window_dphi { {-0.25, 0., 0.},  {0.05, 0., 0.} };
  WindowMatcher window_deta { {100.,0.,0.}, {0.050, 0.0064, 1.1}, {100,0.,0.,}, {0.045, 0.0031, 1.0} };

  bool _print_windows = false;
  void print_windows(bool print=true) { _print_windows = print; }

  void zeroField(const bool flag) { _zero_field = flag; }

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

  TFile *_file = nullptr;
  TNtuple *_tree = nullptr;

  // default values, can be replaced from the macro
  double _phi_search_win = 0.01;
  double _eta_search_win = 0.004;
  double _x_search_win = 0.3;
  double _y_search_win = 0.3;
  double _z_search_win = 0.4;

  bool _use_old_matching = false;  // normally false

  bool _zero_field = false;     // fit straight lines if true

  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeedContainer *_track_map{nullptr};
  TrackSeedContainer *_track_map_silicon{nullptr};
  TrackSeed *_tracklet_tpc{nullptr};
  TrackSeed *_tracklet_si{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  ActsGeometry *_tGeometry{nullptr};
  TrkrClusterCrossingAssoc *_cluster_crossing_map{nullptr};
  int m_event = 0;
  std::map<unsigned int, double> _z_mismatch_map;

  TpcClusterZCrossingCorrection _clusterCrossingCorrection;
  float _crossing_deltaz_max = 10.0;

  //  double _collision_rate = 50e3;  // input rate for phi correction
  //  double _reference_collision_rate = 50e3;  // reference rate for phi correction
  //  double _si_vertex_dzmax = 0.25;  // mm
  double fieldstrength{std::numeric_limits<double>::quiet_NaN()};

  bool _test_windows = false;
  bool _pp_mode = false;
  bool _use_intt_crossing = true;  // should always be true except for testing

  int _n_iteration = 0;
  std::string _track_map_name = "TpcTrackSeedContainer";
  std::string _silicon_track_map_name = "SiliconTrackSeedContainer";
  std::string m_fieldMap = "1.4";
  std::vector<TrkrDefs::cluskey> getTrackletClusterList(TrackSeed* tracklet);
};

#endif  //  PHSILICONTPCTRACKMATCHING_H
