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
  bool _use_legacy_windowing = false;
  void set_use_legacy_windowing (bool set_par=true) { _use_legacy_windowing=set_par; }

  void set_phi_search_window(const double win) { _phi_search_win = win; }
  void set_eta_search_window(const double win) { _eta_search_win = win; }
  void set_x_search_window(const double win) { _x_search_win = win; }
  void set_y_search_window(const double win) { _y_search_win = win; }
  void set_z_search_window(const double win) { _z_search_win = win; }

  float get_phi_search_window() const { return _phi_search_win; }
  float get_eta_search_window() const { return _eta_search_win; }
  float get_x_search_window() const { return _x_search_win; }
  float get_y_search_window() const { return _y_search_win; }
  float get_z_search_window() const { return _z_search_win; }

  // 2024/01/22 update
  // new matching windows are dynamic as functions of Q/pT:
  //   neg Q: a0+b0*exp(c0/pT) < dX < a1+b1*exp(c1/pT)
  //   pos Q: a2+b2*exp(c2/pT) < dX < a3+b3*exp(c3/pT)
  // for x, y:
  //    |dX| < a1       (same for all Q)
  // for phi:
  //     a0 < dphi < a1 (same for all Q)
  // for eta and z, it is:
  //     neg Q: |dX| < a1+b1*exp(c1/pT)
  //     pos Q: |dX| < a3+b3*exp(c3/pT)
  //
  // These are done with this local struct with this logic:
  //    - if a0==100 or a2==100, then do |dX| < a+b*exp(c/pT) for pos and negQ
  //    - if bi (where i=0,1,2,3)==0 don't calculate the exp
  //    - if a1==100, then treat all tracks as positive tracks
  //    - if a3==100, then treat use legacy windowing
  struct WindowMatcher {
    // use legacy values
    bool use_legacy = false;
    double leg_search_win = 1.; // use if use_legacy == true; set in InitRun
    PHSiliconTpcTrackMatching* parent_ptr {nullptr};
    void set_use_legacy(double _leg_search_win, PHSiliconTpcTrackMatching* parent) 
    { 
      use_legacy=true; 
      leg_search_win=_leg_search_win; 
      parent_ptr=parent; 
    }

    using Arr3D = std::array<double,3>;
    Arr3D posR { 100., 0., 0. }; // (above a3,b3,c3), 100 for use _use_legacy_windowing
    Arr3D posL { 100., 0., 0. }; // (above a2,b2,c2), 100 for |dX|
    Arr3D negL { 100., 0., 0. }; // (above a0,b0,c0), 100 for |dX|
    Arr3D negR { 100., 0., 0. }; // (above a1,b1,c1), 100 for treat all tracks pos Q
                                 
    // efficiency flags set during PHSiliconTpcTrackMatching::InitRun()
    double min_pt  = 0.15; // only grow function windows down to 150 MeV
    bool all_pos_Q  = true;
    bool only_fabs  = true;
    bool negL_b0 = true;
    bool negR_b0 = true;
    bool posL_b0 = true;
    bool posR_b0 = true;

    inline double fn_exp(const Arr3D& arr, const bool& b0, double pT) {
      if (b0) { return arr[0]; }
      if (pT<min_pt) pT = min_pt;
      return arr[0]+arr[0]*exp(arr[2]/pT);
    }

    void init_bools();

    bool in_window(bool posQ, double tpc, double si);
    
    WindowMatcher(const Arr3D _posR, const double _min_pt=0.15) 
      : posR{_posR}, min_pt{_min_pt} {};
    WindowMatcher(const Arr3D _posL, const Arr3D _posR, 
        const double _min_pt=0.15) 
      : posR{_posR}, posL{_posL}, min_pt{_min_pt} {};
    WindowMatcher(const Arr3D _posL, const Arr3D _posR, 
         const Arr3D _negL, const Arr3D _negR, const double _min_pt=0.15)
      : posR{_posR}, posL{_posL}, negL{_negL}, negR{_negR}, min_pt{_min_pt} {};

    void set_parameters(const Arr3D _posR, const double _min_pt=0.15) 
    { posR = _posR; min_pt=_min_pt;};
    void set_parameters(const Arr3D _posL, const Arr3D _posR, const double _min_pt=0.15) 
    { posL=_posL; posR=_posR; min_pt=_min_pt;};
    void set_parameters(const Arr3D _posL, const Arr3D _posR, const Arr3D _negL, 
        const Arr3D _negR, const double _min_pt=0.15) 
    { posL=_posL; posR=_posR; negL=_negL; negR=_negR; min_pt=_min_pt;};
  };

  WindowMatcher window_dx   { {5.3, 0., 0.} };
  WindowMatcher window_dy   { {5.2, 0., 0.} };
  WindowMatcher window_dz   { {0., 1.45, 0.49}, {0., 2.6, 0.38} };
  WindowMatcher window_dphi { {-0.25, 0., 0.},  {0.05, 0., 0.} };
  WindowMatcher window_deta { {0.045, 0.0031, 1.0}, {0.050, 0.0064, 1.1} };

  void zeroField(const bool flag) { _zero_field = flag; }
  
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

  TFile *_file = nullptr;
  TNtuple *_tree = nullptr;

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
  std::vector<TrkrDefs::cluskey> getTrackletClusterList(TrackSeed* tracklet);
};

#endif  //  PHSILICONTPCTRACKMATCHING_H
