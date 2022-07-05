#ifndef TRACKBASEHISTORIC_SVTXTRACKV3_H
#define TRACKBASEHISTORIC_SVTXTRACKV3_H

#include "SvtxTrack.h"
#include "SvtxTrackState.h"

#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <utility>  // for pair

class PHObject;

class SvtxTrack_v3: public SvtxTrack
{
 public:
  SvtxTrack_v3();
  
  //* base class copy constructor
  SvtxTrack_v3( const SvtxTrack& );
  
  //* copy constructor
  SvtxTrack_v3(const SvtxTrack_v3& );
  
  //* assignment operator
  SvtxTrack_v3& operator=(const SvtxTrack_v3& track);

  //* destructor
  ~SvtxTrack_v3() override; 

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxTrack_v3(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new SvtxTrack_v3(*this); }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  // copy content from base class
  void CopyFrom( const SvtxTrack& ) override;
  void CopyFrom( SvtxTrack* source ) override
  { CopyFrom( *source ); }

  //
  // basic track information ---------------------------------------------------
  //

  unsigned int get_id() const override { return _track_id; }
  void set_id(unsigned int id) override { _track_id = id; }

  short int get_crossing() const override { return _track_crossing; }
  void set_crossing(short int cross) override { _track_crossing = cross; }

  unsigned int get_vertex_id() const override { return _vertex_id; }
  void set_vertex_id(unsigned int id) override { _vertex_id = id; }

  bool get_positive_charge() const override { return _is_positive_charge; }
  void set_positive_charge(bool ispos) override { _is_positive_charge = ispos; }

  int get_charge() const override { return (get_positive_charge()) ? 1 : -1; }
  void set_charge(int charge) override { (charge > 0) ? set_positive_charge(true) : set_positive_charge(false); }

  float get_chisq() const override { return _chisq; }
  void set_chisq(float chisq) override { _chisq = chisq; }

  unsigned int get_ndf() const override { return _ndf; }
  void set_ndf(int ndf) override { _ndf = ndf; }

  float get_quality() const override { return (_ndf != 0) ? _chisq / _ndf : NAN; }

  float get_dca() const override { return _dca; }
  void set_dca(float dca) override { _dca = dca; }

  float get_dca_error() const override { return _dca_error; }
  void set_dca_error(float dca_error) override { _dca_error = dca_error; }

  float get_dca2d() const override { return _dca2d; }
  void set_dca2d(float dca2d) override { _dca2d = dca2d; }

  float get_dca2d_error() const override { return _dca2d_error; }
  void set_dca2d_error(float error) override { _dca2d_error = error; }

  float get_dca3d_xy() const override { return _dca3d_xy; }
  void set_dca3d_xy(float dcaxy) override { _dca3d_xy = dcaxy; }

  float get_dca3d_xy_error() const override { return _dca3d_xy_error; }
  void set_dca3d_xy_error(float error) override { _dca3d_xy_error = error; }

  float get_dca3d_z() const override { return _dca3d_z; }
  void set_dca3d_z(float dcaz) override { _dca3d_z = dcaz; }

  float get_dca3d_z_error() const override { return _dca3d_z_error; }
  void set_dca3d_z_error(float error) override { _dca3d_z_error = error; }

  float get_x() const override { return _states.find(0.0)->second->get_x(); }
  void set_x(float x) override { _states[0.0]->set_x(x); }

  float get_y() const override { return _states.find(0.0)->second->get_y(); }
  void set_y(float y) override { _states[0.0]->set_y(y); }

  float get_z() const override { return _states.find(0.0)->second->get_z(); }
  void set_z(float z) override { _states[0.0]->set_z(z); }

  float get_pos(unsigned int i) const override { return _states.find(0.0)->second->get_pos(i); }

  float get_px() const override { return _states.find(0.0)->second->get_px(); }
  void set_px(float px) override { _states[0.0]->set_px(px); }

  float get_py() const override { return _states.find(0.0)->second->get_py(); }
  void set_py(float py) override { _states[0.0]->set_py(py); }

  float get_pz() const override { return _states.find(0.0)->second->get_pz(); }
  void set_pz(float pz) override { _states[0.0]->set_pz(pz); }

  float get_mom(unsigned int i) const override { return _states.find(0.0)->second->get_mom(i); }

  float get_p() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2) + pow(get_pz(), 2)); }
  float get_pt() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2)); }
  float get_eta() const override { return asinh(get_pz() / get_pt()); }
  float get_phi() const override { return atan2(get_py(), get_px()); }

  float get_error(int i, int j) const override { return _states.find(0.0)->second->get_error(i, j); }
  void set_error(int i, int j, float value) override { return _states[0.0]->set_error(i, j, value); }

  //
  // state methods -------------------------------------------------------------
  //
  bool empty_states() const override { return _states.empty(); }
  size_t size_states() const override { return _states.size(); }
  size_t count_states(float pathlength) const override { return _states.count(pathlength); }
  void clear_states() override;

  const SvtxTrackState* get_state(float pathlength) const override;
  SvtxTrackState* get_state(float pathlength) override;
  SvtxTrackState* insert_state(const SvtxTrackState* state) override;
  size_t erase_state(float pathlength) override;

  ConstStateIter begin_states() const override { return _states.begin(); }
  ConstStateIter find_state(float pathlength) const override { return _states.find(pathlength); }
  ConstStateIter end_states() const override { return _states.end(); }

  StateIter begin_states() override { return _states.begin(); }
  StateIter find_state(float pathlength) override { return _states.find(pathlength); }
  StateIter end_states() override { return _states.end(); }

  //
  // associated cluster ids methods --------------------------------------------
  //

  // needed by old tracking
  void clear_clusters() override { _cluster_ids.clear(); }
  bool empty_clusters() const override { return _cluster_ids.empty(); }
  size_t size_clusters() const override { return _cluster_ids.size(); }

  void insert_cluster(unsigned int clusterid) override { _cluster_ids.insert(clusterid); }
  size_t erase_cluster(unsigned int clusterid) override { return _cluster_ids.erase(clusterid); }
  ConstClusterIter begin_clusters() const override { return _cluster_ids.begin(); }
  ConstClusterIter find_cluster(unsigned int clusterid) const override { return _cluster_ids.find(clusterid); }
  ConstClusterIter end_clusters() const override { return _cluster_ids.end(); }
  ClusterIter find_cluster(unsigned int clusterid) override { return _cluster_ids.find(clusterid); }
  ClusterIter begin_clusters() override { return _cluster_ids.begin(); }
  ClusterIter end_clusters() override { return _cluster_ids.end(); }

  // needed by new tracking
  void clear_cluster_keys() override { _cluster_keys.clear(); }
  bool empty_cluster_keys() const override { return _cluster_keys.empty(); }
  size_t size_cluster_keys() const override { return _cluster_keys.size(); }

  void insert_cluster_key(TrkrDefs::cluskey clusterid) override { _cluster_keys.insert(clusterid); }
  size_t erase_cluster_key(TrkrDefs::cluskey clusterid) override { return _cluster_keys.erase(clusterid); }
  ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const override { return _cluster_keys.find(clusterid); }
  ConstClusterKeyIter begin_cluster_keys() const override { return _cluster_keys.begin(); }
  ConstClusterKeyIter end_cluster_keys() const override { return _cluster_keys.end(); }
  ClusterKeyIter find_cluster_keys(unsigned int clusterid) override { return _cluster_keys.find(clusterid); }
  ClusterKeyIter begin_cluster_keys() override { return _cluster_keys.begin(); }
  ClusterKeyIter end_cluster_keys() override { return _cluster_keys.end(); }

  //
  // calo projection methods ---------------------------------------------------
  //
  float get_cal_dphi(CAL_LAYER layer) const override;
  void set_cal_dphi(CAL_LAYER layer, float dphi) override { _cal_dphi[layer] = dphi; }

  float get_cal_deta(CAL_LAYER layer) const override;
  void set_cal_deta(CAL_LAYER layer, float deta) override { _cal_deta[layer] = deta; }

  float get_cal_energy_3x3(CAL_LAYER layer) const override;
  void set_cal_energy_3x3(CAL_LAYER layer, float energy_3x3) override { _cal_energy_3x3[layer] = energy_3x3; }

  float get_cal_energy_5x5(CAL_LAYER layer) const override;
  void set_cal_energy_5x5(CAL_LAYER layer, float energy_5x5) override { _cal_energy_5x5[layer] = energy_5x5; }

  unsigned int get_cal_cluster_id(CAL_LAYER layer) const override;
  void set_cal_cluster_id(CAL_LAYER layer, unsigned int id) override { _cal_cluster_id[layer] = id; }

  TrkrDefs::cluskey get_cal_cluster_key(CAL_LAYER layer) const override;
  void set_cal_cluster_key(CAL_LAYER layer, TrkrDefs::cluskey id) override { _cal_cluster_key[layer] = id; }

  float get_cal_cluster_e(CAL_LAYER layer) const override;
  void set_cal_cluster_e(CAL_LAYER layer, float e) override { _cal_cluster_e[layer] = e; }

  // ACTS track information for use by ACTS modules only
  float get_acts_covariance(unsigned int i, unsigned int j) const override {return _acts_trajectory_covariance[i][j]; }
  void set_acts_covariance(unsigned int i, unsigned int j, float value) override { _acts_trajectory_covariance[i][j] = value; }

 private:

  //! acts covariance matrix
  float _acts_trajectory_covariance[6][6] = {};

  // track information
  unsigned int _track_id = UINT_MAX;
  unsigned int _vertex_id = UINT_MAX;
  bool _is_positive_charge = false;
  float _chisq = NAN;
  unsigned int _ndf = 0;
  short int _track_crossing = SHRT_MAX;

  // extended track information (non-primary tracks only)
  float _dca = NAN;
  float _dca_error = NAN;
  float _dca2d = NAN;
  float _dca2d_error = NAN;
  float _dca3d_xy = NAN;
  float _dca3d_xy_error = NAN;
  float _dca3d_z = NAN;
  float _dca3d_z_error = NAN;

  // extended track information (primary tracks only)

  // track state information
  StateMap _states;  //< path length => state object

  // cluster contents
  ClusterSet _cluster_ids;
  ClusterKeySet _cluster_keys;

  // calorimeter matches
  std::map<CAL_LAYER, float> _cal_dphi;
  std::map<CAL_LAYER, float> _cal_deta;
  std::map<CAL_LAYER, float> _cal_energy_3x3;
  std::map<CAL_LAYER, float> _cal_energy_5x5;
  std::map<CAL_LAYER, int> _cal_cluster_id;
  std::map<CAL_LAYER, TrkrDefs::cluskey> _cal_cluster_key;
  std::map<CAL_LAYER, float> _cal_cluster_e;

  ClassDefOverride(SvtxTrack_v3, 2)
};

#endif
