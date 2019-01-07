#ifndef __SVTXVERTEX_V1_H__
#define __SVTXVERTEX_V1_H__

#include "SvtxVertex.h"

#include <phool/PHObject.h>
#include <iostream>
#include <set>
#include <vector>

class SvtxVertex_v1 : public SvtxVertex
{
 public:
  SvtxVertex_v1();
  virtual ~SvtxVertex_v1() {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const;
  void Reset() { *this = SvtxVertex_v1(); }
  int isValid() const;
  SvtxVertex* Clone() const { return new SvtxVertex_v1(*this); }

  // vertex info

  unsigned int get_id() const { return _id; }
  void set_id(unsigned int id) { _id = id; }

  float get_t0() const { return _t0; }
  void set_t0(float t0) { _t0 = t0; }

  float get_x() const { return _pos[0]; }
  void set_x(float x) { _pos[0] = x; }

  float get_y() const { return _pos[1]; }
  void set_y(float y) { _pos[1] = y; }

  float get_z() const { return _pos[2]; }
  void set_z(float z) { _pos[2] = z; }

  float get_chisq() const { return _chisq; }
  void set_chisq(float chisq) { _chisq = chisq; }

  unsigned int get_ndof() const { return _ndof; }
  void set_ndof(unsigned int ndof) { _ndof = ndof; }

  float get_position(unsigned int coor) const { return _pos[coor]; }
  void set_position(unsigned int coor, float xi) { _pos[coor] = xi; }

  float get_error(unsigned int i, unsigned int j) const;        //< get vertex error covar
  void set_error(unsigned int i, unsigned int j, float value);  //< set vertex error covar

  //
  // associated track ids methods
  //
  void clear_tracks() { _track_ids.clear(); }
  bool empty_tracks() { return _track_ids.empty(); }
  size_t size_tracks() const { return _track_ids.size(); }
  void insert_track(unsigned int trackid) { _track_ids.insert(trackid); }
  size_t erase_track(unsigned int trackid) { return _track_ids.erase(trackid); }
  ConstTrackIter begin_tracks() const { return _track_ids.begin(); }
  ConstTrackIter find_track(unsigned int trackid) const { return _track_ids.find(trackid); }
  ConstTrackIter end_tracks() const { return _track_ids.end(); }
  TrackIter begin_tracks() { return _track_ids.begin(); }
  TrackIter find_track(unsigned int trackid) { return _track_ids.find(trackid); }
  TrackIter end_tracks() { return _track_ids.end(); }

 private:
  unsigned int covar_index(unsigned int i, unsigned int j) const;

  unsigned int _id;                   //< unique identifier within container
  float _t0;                          //< collision time
  float _pos[3];                      //< collision position x,y,z
  float _chisq;                       //< vertex fit chisq
  unsigned int _ndof;                 //< degrees of freedom
  float _err[6];                      //< error covariance matrix (packed storage) (+/- cm^2)
  std::set<unsigned int> _track_ids;  //< list of track ids

  ClassDef(SvtxVertex_v1, 1);
};

#endif
