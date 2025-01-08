#ifndef TRACKBASEHISTORIC_SVTXVERTEX_H
#define TRACKBASEHISTORIC_SVTXVERTEX_H

#include "Vertex.h"

#include <climits>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

class SvtxVertex : public Vertex
{
 public:
  ~SvtxVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxVertex base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  // vertex info

  virtual unsigned int get_id() const override { return UINT_MAX; }
  virtual void set_id(unsigned int) override {}

  virtual float get_t0() const override { return get_t(); }
  virtual void set_t0(float t0) override { set_t(t0); }

  virtual float get_t() const override { return NAN; }
  virtual void set_t(float) override {}

  virtual float get_x() const override { return NAN; }
  virtual void set_x(float) override {}

  virtual float get_y() const override { return NAN; }
  virtual void set_y(float) override {}

  virtual float get_z() const override { return NAN; }
  virtual void set_z(float) override {}

  virtual float get_chisq() const override { return NAN; }
  virtual void set_chisq(float) override {}

  virtual unsigned int get_ndof() const override { return UINT_MAX; }
  virtual void set_ndof(unsigned int) override {}

  virtual float get_position(unsigned int) const override { return NAN; }
  virtual void set_position(unsigned int, float) override {}

  virtual float get_error(unsigned int, unsigned int) const override { return NAN; }
  virtual void set_error(unsigned int, unsigned int, float) override {}

  //
  // associated track ids methods
  //
  virtual void clear_tracks() override {}
  virtual bool empty_tracks() override { return true; }
  virtual size_t size_tracks() const override { return 0; }
  virtual void insert_track(unsigned int /*trackid*/) override {}
  virtual size_t erase_track(unsigned int /*trackid*/) override { return 0; }
  virtual ConstTrackIter begin_tracks() const override;
  virtual ConstTrackIter find_track(unsigned int trackid) const override;
  virtual ConstTrackIter end_tracks() const override;
  virtual TrackIter begin_tracks() override;
  virtual TrackIter find_track(unsigned int trackid) override;
  virtual TrackIter end_tracks() override;

 protected:
  ClassDefOverride(SvtxVertex, 1);
};

#endif
