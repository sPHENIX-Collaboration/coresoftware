#ifndef TRACKBASEHISTORIC_SVTXTRACKV4_H
#define TRACKBASEHISTORIC_SVTXTRACKV4_H

#include "SvtxTrack.h"
#include "SvtxTrackState.h"
#include "TrackSeed.h"

#include <trackbase/TrkrDefs.h>

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <utility>  // for pair

class PHObject;

class SvtxTrack_v4: public SvtxTrack
{
 public:
  SvtxTrack_v4();
  
  //* base class copy constructor
  SvtxTrack_v4( const SvtxTrack& );
  
  //* copy constructor
  SvtxTrack_v4(const SvtxTrack_v4& );
  
  //* assignment operator
  SvtxTrack_v4& operator=(const SvtxTrack_v4& track);

  //* destructor
  ~SvtxTrack_v4() override; 

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxTrack_v4(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new SvtxTrack_v4(*this); }

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

  TrackSeed* get_tpc_seed() const override { return _tpc_seed; }
  void set_tpc_seed(TrackSeed* seed) override { _tpc_seed = seed; }

  TrackSeed* get_silicon_seed() const override { return _silicon_seed; }
  void set_silicon_seed(TrackSeed* seed) override {_silicon_seed = seed; }

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
 
 private:

  // track information
  TrackSeed* _tpc_seed = nullptr;
  TrackSeed* _silicon_seed = nullptr;
  unsigned int _track_id = UINT_MAX;
  unsigned int _vertex_id = UINT_MAX;
  bool _is_positive_charge = false;
  float _chisq = NAN;
  unsigned int _ndf = 0;
  short int _track_crossing = SHRT_MAX;

  // track state information
  StateMap _states;  //< path length => state object

  ClassDefOverride(SvtxTrack_v4, 4)
};

#endif
