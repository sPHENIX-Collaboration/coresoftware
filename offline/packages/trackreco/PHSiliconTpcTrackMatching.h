// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONTPCTRACKMATCHING_H
#define PHSILICONTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <map>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class TF1;
class TpcSeedTrackMap;
class AssocInfoContainer;
class TrkrClusterContainer;

class PHSiliconTpcTrackMatching : public SubsysReco
{
 public:

  PHSiliconTpcTrackMatching(const std::string &name = "PHSiliconTpcTrackMatching");

  ~PHSiliconTpcTrackMatching() override;

  void set_track_map_name_silicon(const std::string &map_name) { _track_map_name_silicon = map_name; }
  void set_phi_search_window(const double win){_phi_search_win = win;}
  void set_eta_search_window(const double win){_eta_search_win = win;}
  void set_x_search_window(const double win){_x_search_win = win;}
  void set_y_search_window(const double win){_y_search_win = win;}
  void set_z_search_window(const double win){_z_search_win = win;}
  void set_search_par_values(const double p0, const double p1, const double p2){_par0 = p0; _par1 = p1; _par2 = p2; }
  void set_seeder(const bool is_ca_seeder){_is_ca_seeder = is_ca_seeder;}

  void set_field_dir(const double rescale)
  {
    _fieldDir = -1;
    if(rescale > 0)
      _fieldDir = 1;     
  }
  void set_field(const std::string &field) { _field = field;}

  void set_test_windows_printout(const bool test){_test_windows = test ;}
  void set_sc_calib_mode(const bool flag){_sc_calib_flag = flag;}
  void set_collision_rate(const double rate){_collision_rate = rate;}

  int InitRun(PHCompositeNode* topNode) override;

  int process_event(PHCompositeNode*) override;

  int End(PHCompositeNode*) override;

  void set_silicon_track_map_name(const std::string &map_name) { _silicon_track_map_name = map_name; }
  void set_tpcseed_track_map_name(const std::string &map_name) { _tpcseed_track_map_name = map_name; }
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void SetIteration(int iter){_n_iteration = iter;}
 private:

  int GetNodes(PHCompositeNode* topNode);

  void copySiliconClustersToCorrectedMap( );

  std::string _track_map_name_silicon;

  // default values, can be replaced from the macro
  double _phi_search_win = 0.01;
  double _eta_search_win = 0.004;
  double _x_search_win = 0.3;
  double _y_search_win = 0.3;
  double _z_search_win = 0.4;
  
  AssocInfoContainer *_assoc_container{nullptr};
  SvtxTrackMap *_track_map{nullptr};
  SvtxTrackMap *_track_map_silicon{nullptr};
  SvtxTrack *_tracklet_tpc{nullptr};
  SvtxTrack *_tracklet_si{nullptr};
  TrkrClusterContainer *_cluster_map{nullptr};
  TrkrClusterContainer *_corrected_cluster_map{nullptr};

  TpcSeedTrackMap *_seed_track_map{nullptr};
  //std::multimap<unsigned int, unsigned int> _seed_track_map;
 
  // correction function for PHTpcTracker track phi bias
  TF1 *fdphi{nullptr};
  // default values, can be replaced from the macro
  double _par0 =  -0.000650;
  double _par1 =  0.13373;
  double _par2 =  0.98298;

  // correction factor for TPC tracklet phi offset due to space charge in TPC
  TF1 *fscdphi{nullptr};
  double _parsc0 =  0.0366293;
  double _parsc1 =   -0.0133073;

  double _collision_rate = 50e3;  // input rate for phi correction
  double _reference_collision_rate = 50e3;  // reference rate for phi correction

  bool _is_ca_seeder = true;
  bool _sc_calib_flag = false;
  bool _test_windows = false;

  std::string _field;
  int _fieldDir = -1;

  int _n_iteration = 0;
  std::string _track_map_name = "SvtxTrackMap";
  std::string _tpcseed_track_map_name = "TpcSeedTrackMap";
  std::string _silicon_track_map_name = "SvtxSiliconTrackMap";
};

#endif // PHTRUTHSILICONASSOCIATION_H
