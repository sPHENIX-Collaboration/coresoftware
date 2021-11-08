#ifndef TRACKBASEHISTORIC_SVTXVERTEX_H
#define TRACKBASEHISTORIC_SVTXVERTEX_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <climits>
#include <set>
#include <vector>

class SvtxVertex : public PHObject
{
 public:
  typedef std::set<unsigned int> TrackSet;
  typedef std::set<unsigned int>::const_iterator ConstTrackIter;
  typedef std::set<unsigned int>::iterator TrackIter;

  ~SvtxVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxVertex base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  // vertex info

  virtual unsigned int get_id() const { return UINT_MAX; }
  virtual void set_id(unsigned int) {}

  virtual float get_t0() const { return NAN; }
  virtual void set_t0(float) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_chisq() const { return NAN; }
  virtual void set_chisq(float) {}

  virtual unsigned int get_ndof() const { return UINT_MAX; }
  virtual void set_ndof(unsigned int) {}

  virtual float get_position(unsigned int) const { return NAN; }
  virtual void set_position(unsigned int /*coor*/, float /*xi*/) {}

  virtual float get_error(unsigned int /*i*/, unsigned int /*j*/) const { return NAN; }
  virtual void set_error(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}

  //
  // associated track ids methods
  //
  virtual void clear_tracks() {}
  virtual bool empty_tracks() { return true; }
  virtual size_t size_tracks() const { return 0; }
  virtual void insert_track(unsigned int /*trackid*/) {}
  virtual size_t erase_track(unsigned int /*trackid*/) { return 0; }
  virtual ConstTrackIter begin_tracks() const;
  virtual ConstTrackIter find_track(unsigned int trackid) const;
  virtual ConstTrackIter end_tracks() const;
  virtual TrackIter begin_tracks();
  virtual TrackIter find_track(unsigned int trackid);
  virtual TrackIter end_tracks();

 protected:
  SvtxVertex() {}

  ClassDefOverride(SvtxVertex, 1);
};

#endif
