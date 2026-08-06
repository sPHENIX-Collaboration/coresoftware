#ifndef TRACKBASEHISTORIC_SVTXTRACKV5_H
#define TRACKBASEHISTORIC_SVTXTRACKV5_H

#include "SvtxTrack.h"
#include "SvtxTrackState.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>

class PHObject;

class SvtxTrack_v5 : public SvtxTrack
{
 public:
  SvtxTrack_v5();

  //* base class copy constructor
  SvtxTrack_v5(const SvtxTrack&);

  //* copy constructor
  SvtxTrack_v5(const SvtxTrack_v5&);

  //* assignment operator
  SvtxTrack_v5& operator=(const SvtxTrack_v5& source);

  //* destructor
  ~SvtxTrack_v5() override;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxTrack_v5(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new SvtxTrack_v5(*this); }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;
  // copy content from base class
  void CopyFrom(const SvtxTrack&) override;
  void CopyFrom(SvtxTrack* source) override
  {
    CopyFrom(*source);
  }

  //
  // basic track information ---------------------------------------------------
  //

  unsigned int get_id() const override { return _track_id; }
  void set_id(unsigned int id) override { _track_id = id; }

  short int get_crossing() const override { return _track_crossing; }
  void set_crossing(short int crossing) override { _track_crossing = crossing; }

  unsigned int get_vertex_id() const override { return _vertex_id; }
  void set_vertex_id(unsigned int id) override { _vertex_id = id; }

  bool get_positive_charge() const override { return m_charge > 0; }
  void set_positive_charge(bool ispos) override { m_charge = ispos ? 1 : -1; }

  int get_charge() const override { return m_charge; }
  void set_charge(int charge) override { m_charge = (charge > 0) ? 1 : -1; }

  float get_chisq() const override { return _chisq; }
  void set_chisq(float chisq) override { _chisq = chisq; }

  unsigned int get_ndf() const override { return _ndf; }
  void set_ndf(int ndf) override;

  float get_quality() const override { return (_ndf != 0) ? _chisq / _ndf : NAN; }

  float get_x() const override;
  void set_x(float x) override;

  float get_y() const override;
  void set_y(float y) override;

  float get_z() const override;
  void set_z(float z) override;

  float get_pos(unsigned int i) const override;

  float get_px() const override;
  void set_px(float px) override;

  float get_py() const override;
  void set_py(float py) override;

  float get_pz() const override;
  void set_pz(float pz) override;

  float get_mom(unsigned int i) const override;

  float get_p() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2) + pow(get_pz(), 2)); }
  float get_pt() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2)); }
  float get_eta() const override { return asinh(get_pz() / get_pt()); }
  float get_phi() const override { return atan2(get_py(), get_px()); }

  float get_error(int i, int j) const override;
  void set_error(int i, int j, float value) override;

  //
  // cluster count information -------------------------------------------------
  //

  unsigned int get_nmvtx_clusters() const { return _nmvtx_clusters; }
  void set_nmvtx_clusters(unsigned int nclusters) { _nmvtx_clusters = compress_cluster_count(nclusters); }

  unsigned int get_nintt_clusters() const { return _nintt_clusters; }
  void set_nintt_clusters(unsigned int nclusters) { _nintt_clusters = compress_cluster_count(nclusters); }

  unsigned int get_ntpc_clusters() const { return _ntpc_clusters; }
  void set_ntpc_clusters(unsigned int nclusters) { _ntpc_clusters = compress_cluster_count(nclusters); }

  unsigned int get_ntpot_clusters() const { return _ntpot_clusters; }
  void set_ntpot_clusters(unsigned int nclusters) { _ntpot_clusters = compress_cluster_count(nclusters); }

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

 private:
  const SvtxTrackState* get_pca_state() const;
  SvtxTrackState* get_pca_state();

  // track state information
  StateMap _states;  //< path length => state object. each state has 30 floats, a string, and a cluskey

  // track information
  float _chisq = std::numeric_limits<float>::quiet_NaN(); //4byte
  unsigned int _track_id = std::numeric_limits<unsigned int>::max(); //4byte
  unsigned int _vertex_id = std::numeric_limits<unsigned int>::max(); //4byte
  short int _track_crossing = std::numeric_limits<short int>::max(); //2byte
  signed char m_charge = 0; //1byte
  unsigned char _ndf = 0; //1byte
  unsigned char _nmvtx_clusters = 0; //1byte
  unsigned char _nintt_clusters = 0; //1byte
  unsigned char _ntpc_clusters = 0; //1byte
  unsigned char _ntpot_clusters = 0; //1byte

  ClassDefOverride(SvtxTrack_v5, 1)
};

#endif
