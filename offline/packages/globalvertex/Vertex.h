#ifndef VERTEX_H
#define VERTEX_H

#include <phool/PHObject.h>
#include <climits>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

class Vertex : public PHObject
{
 public:
  typedef std::set<unsigned int> TrackSet;
  typedef std::set<unsigned int>::const_iterator ConstTrackIter;
  typedef std::set<unsigned int>::iterator TrackIter;

  ~Vertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "Vertex base class" << std::endl; }
  PHObject* CloneMe() const override { return nullptr; }
  int isValid() const override { return 0; }

  // vertex info

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) {}

  virtual float get_t() const { return NAN; }
  virtual void set_t(float) {}

  // Interface functions to maintain backwards compatibility with svtxvertex_v1
  virtual float get_t0() const { return get_t(); }
  virtual void set_t0(float t0) { set_t(t0); }

  virtual float get_t_err() const { return NAN; }
  virtual void set_t_err(float) {}

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

  // beam crossing methods
  virtual unsigned int get_beam_crossing() const { return UINT_MAX; }
  virtual void set_beam_crossing(unsigned int) {}

  // bbcvertex methods
  virtual void set_bbc_ns(int, int, float, float) {}
  virtual int get_bbc_npmt(int) const { return std::numeric_limits<int>::max(); }
  virtual float get_z_err() const { return NAN; }
  virtual void set_z_err(float) {}

  virtual float get_bbc_q(int) const { return NAN; }
  virtual float get_bbc_t(int) const { return NAN; }

  // svtxvertex methods
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

 private:
  ClassDefOverride(Vertex, 1);
};

#endif
