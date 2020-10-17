// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONTPCTRACKMATCHING_H
#define PHSILICONTPCTRACKMATCHING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <set>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class PHG4Particle;
class AssocInfoContainer;
class TF1;

#include <trackreco/PHTrackPropagating.h>

class PHSiliconTpcTrackMatching : public PHTrackPropagating
{
 public:

  PHSiliconTpcTrackMatching(const std::string &name = "PHSiliconTpcTrackMatching");

  virtual ~PHSiliconTpcTrackMatching();

  void set_track_map_name_silicon(const std::string &map_name) { _track_map_name_silicon = map_name; }
  void set_phi_search_window(const double win){_phi_search_win = win;}
  void set_eta_search_window(const double win){_eta_search_win = win;}
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

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  std::string _track_map_name_silicon;

  // default values, can be replaced from the macro
  double _phi_search_win = 0.02;
  double _eta_search_win = 0.015;
  
  SvtxTrackMap *_track_map_silicon{nullptr};
  SvtxTrack *_tracklet_tpc{nullptr};
  SvtxTrack *_tracklet_si{nullptr};

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

  bool _is_ca_seeder = false;
  bool _sc_calib_flag = false;
  bool _test_windows = false;

  std::string _field;
  int _fieldDir = -1;

};

#endif // PHTRUTHSILICONASSOCIATION_H
