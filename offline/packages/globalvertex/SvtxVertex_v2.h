// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_SVTXVERTEXV2_H
#define GLOBALVERTEX_SVTXVERTEXV2_H

#include "SvtxVertex.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <set>

class PHObject;

class SvtxVertex_v2 : public SvtxVertex
{
 public:
  SvtxVertex_v2();
  ~SvtxVertex_v2() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxVertex_v2(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new SvtxVertex_v2(*this); }

  // vertex info

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  float get_t() const override { return _t0; }
  void set_t(float t0) override { _t0 = t0; }

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

  short int get_beam_crossing() const override 
  { 
    return rollover_from_unsignedint(_beamcrossing);
  }
  void set_beam_crossing(short int cross) override 
  { 
    if (cross == short_int_max)
    {
      _beamcrossing = std::numeric_limits<unsigned int>::max();
      return;
    }

    const short int cross_ro = rollover_short(cross);
    _beamcrossing = static_cast<unsigned int>(cross_ro);
  }

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
  static constexpr short int short_int_max = std::numeric_limits<short int>::max(); // 32767
  // for unsigned int to short int conversion (rollover)
  static short int rollover_short(short int cross)
  {
    if (cross == short_int_max) return short_int_max;
    if (cross >= 0) return cross;

    const int cross_ro = static_cast<int>(short_int_max) + static_cast<int>(cross); // cross negative
    return static_cast<short int>(cross_ro);
  }

  static short int rollover_from_unsignedint(unsigned int cross)
  {
    // if unsigned int max, return short int max
    if (cross == std::numeric_limits<unsigned int>::max())
    {
      return short_int_max;
    }

    // Common case: [0, 32767]
    if (cross <= static_cast<unsigned int>(short_int_max))
    {
      return static_cast<short int>(cross);
    }

    const short int cross_ro = static_cast<short int>(static_cast<unsigned short>(cross));
    if (cross_ro >= 0) return cross_ro;

    return rollover_short(cross_ro);
  }

  unsigned int covar_index(unsigned int i, unsigned int j) const;

  unsigned int _id{std::numeric_limits<unsigned int>::max()};    //< unique identifier within container
  float _t0{std::numeric_limits<float>::quiet_NaN()};            //< collision time
  float _pos[3]{};                                               //< collision position x,y,z
  float _chisq{std::numeric_limits<float>::quiet_NaN()};         //< vertex fit chisq
  unsigned int _ndof{std::numeric_limits<unsigned int>::max()};  //< degrees of freedom
  float _err[6]{};                                               //< error covariance matrix (packed storage) (+/- cm^2)
  std::set<unsigned int> _track_ids;                             //< list of track ids
  unsigned int _beamcrossing{std::numeric_limits<unsigned int>::max()};

  ClassDefOverride(SvtxVertex_v2, 2);
};

#endif
