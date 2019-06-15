#ifndef __SVTXVERTEX_H__
#define __SVTXVERTEX_H__

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <limits.h>
#include <set>
#include <vector>

class SvtxVertex : public PHObject
{
 public:
  typedef std::set<unsigned int> TrackSet;
  typedef std::set<unsigned int>::const_iterator ConstTrackIter;
  typedef std::set<unsigned int>::iterator TrackIter;

  virtual ~SvtxVertex() {}

  // PHObject virtual overloads

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "SvtxVertex base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual SvtxVertex* clone() const { return NULL; }

  // vertex info

  virtual unsigned int get_id() const { return UINT_MAX; }
  virtual void set_id(unsigned int id) {}

  virtual float get_t0() const { return NAN; }
  virtual void set_t0(float t0) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float x) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float y) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float z) {}

  virtual float get_chisq() const { return NAN; }
  virtual void set_chisq(float chisq) {}

  virtual unsigned int get_ndof() const { return UINT_MAX; }
  virtual void set_ndof(unsigned int ndof) {}

  virtual float get_position(unsigned int coor) const { return NAN; }
  virtual void set_position(unsigned int coor, float xi) {}

  virtual float get_error(unsigned int i, unsigned int j) const { return NAN; }
  virtual void set_error(unsigned int i, unsigned int j, float value) {}

  //
  // associated track ids methods
  //
  virtual void clear_tracks() {}
  virtual bool empty_tracks() { return true; }
  virtual size_t size_tracks() const { return 0; }
  virtual void insert_track(unsigned int trackid) {}
  virtual size_t erase_track(unsigned int trackid) { return 0; }
  virtual ConstTrackIter begin_tracks() const { return TrackSet().end(); }
  virtual ConstTrackIter find_track(unsigned int trackid) const { return TrackSet().end(); }
  virtual ConstTrackIter end_tracks() const { return TrackSet().end(); }
  virtual TrackIter begin_tracks() { return TrackSet().end(); }
  virtual TrackIter find_track(unsigned int trackid) { return TrackSet().end(); }
  virtual TrackIter end_tracks() { return TrackSet().end(); }

 protected:
  SvtxVertex() {}

  ClassDef(SvtxVertex, 1);
};

#endif
