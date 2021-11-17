#ifndef TRACKBASEHISTORIC_SVTXVERTEXV1_H
#define TRACKBASEHISTORIC_SVTXVERTEXV1_H

#include "SvtxVertex.h"

#include <cstddef>      // for size_t
#include <iostream>
#include <set>

class PHObject;

class SvtxVertex_v1 : public SvtxVertex
{
 public:
  SvtxVertex_v1();
  ~SvtxVertex_v1() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxVertex_v1(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new SvtxVertex_v1(*this); }

  // vertex info

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  float get_t0() const override { return _t0; }
  void set_t0(float t0) override { _t0 = t0; }

  float get_x() const override { return _pos[0]; }
  void set_x(float x) override { _pos[0] = x; }

  float get_y() const override { return _pos[1]; }
  void set_y(float y) override { _pos[1] = y; }

  float get_z() const override { return _pos[2]; }
  void set_z(float z) override { _pos[2] = z; }

  float get_chisq() const override { return _chisq; }
  void set_chisq(float chisq) override { _chisq = chisq; }

  unsigned int get_ndof() const override { return _ndof; }
  void set_ndof(unsigned int ndof) override { _ndof = ndof; }

  float get_position(unsigned int coor) const override { return _pos[coor]; }
  void set_position(unsigned int coor, float xi) override { _pos[coor] = xi; }

  float get_error(unsigned int i, unsigned int j) const override;        //< get vertex error covar
  void set_error(unsigned int i, unsigned int j, float value) override;  //< set vertex error covar

  //
  // associated track ids methods
  //
  void clear_tracks() override { _track_ids.clear(); }
  bool empty_tracks() override { return _track_ids.empty(); }
  size_t size_tracks() const override { return _track_ids.size(); }
  void insert_track(unsigned int trackid) override { _track_ids.insert(trackid); }
  size_t erase_track(unsigned int trackid) override { return _track_ids.erase(trackid); }
  ConstTrackIter begin_tracks() const override { return _track_ids.begin(); }
  ConstTrackIter find_track(unsigned int trackid) const override { return _track_ids.find(trackid); }
  ConstTrackIter end_tracks() const override { return _track_ids.end(); }
  TrackIter begin_tracks() override { return _track_ids.begin(); }
  TrackIter find_track(unsigned int trackid) override { return _track_ids.find(trackid); }
  TrackIter end_tracks() override { return _track_ids.end(); }

 private:
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  unsigned int _id;                   //< unique identifier within container
  float _t0;                          //< collision time
  float _pos[3];                      //< collision position x,y,z
  float _chisq;                       //< vertex fit chisq
  unsigned int _ndof;                 //< degrees of freedom
  float _err[6];                      //< error covariance matrix (packed storage) (+/- cm^2)
  std::set<unsigned int> _track_ids;  //< list of track ids

  ClassDefOverride(SvtxVertex_v1, 1);
};

#endif
