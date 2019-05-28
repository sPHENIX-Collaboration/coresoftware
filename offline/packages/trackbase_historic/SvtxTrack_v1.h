#ifndef TRACKBASEHISTORIC_SVTXTRACKV1_H
#define TRACKBASEHISTORIC_SVTXTRACKV1_H

#include "SvtxTrack.h"
#include "SvtxTrackState.h"

#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <cstddef>              // for size_t
#include <iostream>
#include <map>
#include <utility>               // for pair

class SvtxTrack_v1 : public SvtxTrack
{
 public:
  SvtxTrack_v1();
  SvtxTrack_v1(const SvtxTrack_v1& track);
  SvtxTrack_v1& operator=(const SvtxTrack_v1& track);
  virtual ~SvtxTrack_v1();

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const;
  void Reset() { *this = SvtxTrack_v1(); }
  int isValid() const;
  SvtxTrack* Clone() const { return new SvtxTrack_v1(*this); }

  //
  // basic track information ---------------------------------------------------
  //

  unsigned int get_id() const { return _track_id; }
  void set_id(unsigned int id) { _track_id = id; }

  bool get_positive_charge() const { return _is_positive_charge; }
  void set_positive_charge(bool ispos) { _is_positive_charge = ispos; }

  int get_charge() const { return (get_positive_charge()) ? 1 : -1; }
  void set_charge(int charge) { (charge > 0) ? set_positive_charge(true) : set_positive_charge(false); }

  float get_chisq() const { return _chisq; }
  void set_chisq(float chisq) { _chisq = chisq; }

  unsigned int get_ndf() const { return _ndf; }
  void set_ndf(int ndf) { _ndf = ndf; }

  float get_quality() const { return (_ndf != 0) ? _chisq / _ndf : NAN; }

  float get_dca() const { return _dca; }
  void set_dca(float dca) { _dca = dca; }

  float get_dca_error() const { return _dca_error; }
  void set_dca_error(float dca_error) { _dca_error = dca_error; }

  float get_dca2d() const { return _dca2d; }
  void set_dca2d(float dca2d) { _dca2d = dca2d; }

  float get_dca2d_error() const { return _dca2d_error; }
  void set_dca2d_error(float error) { _dca2d_error = error; }

  float get_dca3d_xy() const { return _dca3d_xy; }
  void set_dca3d_xy(float dcaxy) { _dca3d_xy = dcaxy; }

  float get_dca3d_xy_error() const { return _dca3d_xy_error; }
  void set_dca3d_xy_error(float error) { _dca3d_xy_error = error; }

  float get_dca3d_z() const { return _dca3d_z; }
  void set_dca3d_z(float dcaz) { _dca3d_z = dcaz; }

  float get_dca3d_z_error() const { return _dca3d_z_error; }
  void set_dca3d_z_error(float error) { _dca3d_z_error = error; }

  float get_x() const { return _states.find(0.0)->second->get_x(); }
  void set_x(float x) { _states[0.0]->set_x(x); }

  float get_y() const { return _states.find(0.0)->second->get_y(); }
  void set_y(float y) { _states[0.0]->set_y(y); }

  float get_z() const { return _states.find(0.0)->second->get_z(); }
  void set_z(float z) { _states[0.0]->set_z(z); }

  float get_pos(unsigned int i) const { return _states.find(0.0)->second->get_pos(i); }

  float get_px() const { return _states.find(0.0)->second->get_px(); }
  void set_px(float px) { _states[0.0]->set_px(px); }

  float get_py() const { return _states.find(0.0)->second->get_py(); }
  void set_py(float py) { _states[0.0]->set_py(py); }

  float get_pz() const { return _states.find(0.0)->second->get_pz(); }
  void set_pz(float pz) { _states[0.0]->set_pz(pz); }

  float get_mom(unsigned int i) const { return _states.find(0.0)->second->get_mom(i); }

  float get_p() const { return sqrt(pow(get_px(), 2) + pow(get_py(), 2) + pow(get_pz(), 2)); }
  float get_pt() const { return sqrt(pow(get_px(), 2) + pow(get_py(), 2)); }
  float get_eta() const { return asinh(get_pz() / get_pt()); }
  float get_phi() const { return atan2(get_py(), get_px()); }

  float get_error(int i, int j) const { return _states.find(0.0)->second->get_error(i, j); }
  void set_error(int i, int j, float value) { return _states[0.0]->set_error(i, j, value); }

  //
  // state methods -------------------------------------------------------------
  //
  bool empty_states() const { return _states.empty(); }
  size_t size_states() const { return _states.size(); }
  size_t count_states(float pathlength) const { return _states.count(pathlength); }
  void clear_states();

  const SvtxTrackState* get_state(float pathlength) const;
  SvtxTrackState* get_state(float pathlength);
  SvtxTrackState* insert_state(const SvtxTrackState* state);
  size_t erase_state(float pathlength);

  ConstStateIter begin_states() const { return _states.begin(); }
  ConstStateIter find_state(float pathlength) const { return _states.find(pathlength); }
  ConstStateIter end_states() const { return _states.end(); }

  StateIter begin_states() { return _states.begin(); }
  StateIter find_state(float pathlength) { return _states.find(pathlength); }
  StateIter end_states() { return _states.end(); }

  //
  // associated cluster ids methods --------------------------------------------
  //

  // needed by old tracking
  void clear_clusters() { _cluster_ids.clear(); }
  bool empty_clusters() const { return _cluster_ids.empty(); }
  size_t size_clusters() const { return _cluster_ids.size(); }

  void insert_cluster(unsigned int clusterid) { _cluster_ids.insert(clusterid); }
  size_t erase_cluster(unsigned int clusterid) { return _cluster_ids.erase(clusterid); }
  ConstClusterIter begin_clusters() const { return _cluster_ids.begin(); }
  ConstClusterIter find_cluster(unsigned int clusterid) const { return _cluster_ids.find(clusterid); }
  ConstClusterIter end_clusters() const { return _cluster_ids.end(); }
  ClusterIter find_cluster(unsigned int clusterid) { return _cluster_ids.find(clusterid); }
  ClusterIter begin_clusters() { return _cluster_ids.begin(); }
  ClusterIter end_clusters() { return _cluster_ids.end(); }

  // needed by new tracking
  void clear_cluster_keys() { _cluster_keys.clear(); }
  bool empty_cluster_keys() const { return _cluster_keys.empty(); }
  size_t size_cluster_keys() const { return _cluster_keys.size(); }

  void insert_cluster_key(TrkrDefs::cluskey clusterid) { _cluster_keys.insert(clusterid); }
  size_t erase_cluster_key(TrkrDefs::cluskey clusterid) { return _cluster_keys.erase(clusterid); }
  ConstClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) const { return _cluster_keys.find(clusterid); }
  ConstClusterKeyIter begin_cluster_keys() const { return _cluster_keys.begin(); }
  ConstClusterKeyIter end_cluster_keys() const { return _cluster_keys.end(); }
  ClusterKeyIter find_cluster_key(TrkrDefs::cluskey clusterid) { return _cluster_keys.find(clusterid); }
  ClusterKeyIter begin_cluster_keys() { return _cluster_keys.begin(); }
  ClusterKeyIter end_cluster_keys()  { return _cluster_keys.end(); }

  //
  // calo projection methods ---------------------------------------------------
  //
  float get_cal_dphi(CAL_LAYER layer) const;
  void set_cal_dphi(CAL_LAYER layer, float dphi) { _cal_dphi[layer] = dphi; }

  float get_cal_deta(CAL_LAYER layer) const;
  void set_cal_deta(CAL_LAYER layer, float deta) { _cal_deta[layer] = deta; }

  float get_cal_energy_3x3(CAL_LAYER layer) const;
  void set_cal_energy_3x3(CAL_LAYER layer, float energy_3x3) { _cal_energy_3x3[layer] = energy_3x3; }

  float get_cal_energy_5x5(CAL_LAYER layer) const;
  void set_cal_energy_5x5(CAL_LAYER layer, float energy_5x5) { _cal_energy_5x5[layer] = energy_5x5; }

  unsigned int get_cal_cluster_id(CAL_LAYER layer) const;
  void set_cal_cluster_id(CAL_LAYER layer, unsigned int id) { _cal_cluster_id[layer] = id; }

  TrkrDefs::cluskey get_cal_cluster_key(CAL_LAYER layer) const;
  void set_cal_cluster_key(CAL_LAYER layer, TrkrDefs::cluskey id) { _cal_cluster_key[layer] = id; }

  float get_cal_cluster_e(CAL_LAYER layer) const;
  void set_cal_cluster_e(CAL_LAYER layer, float e) { _cal_cluster_e[layer] = e; }

 private:
  // track information
  unsigned int _track_id;
  bool _is_positive_charge;
  float _chisq;
  unsigned int _ndf;

  // extended track information (non-primary tracks only)
  float _dca;
  float _dca_error;
  float _dca2d;
  float _dca2d_error;
  float _dca3d_xy;
  float _dca3d_xy_error;
  float _dca3d_z;
  float _dca3d_z_error;

  // extended track information (primary tracks only)
  // unsigned int _vertex_id;

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

  ClassDef(SvtxTrack_v1, 1)
};

#endif
