#ifndef __TRACKERLUSTER_H__
#define __TRACKERLUSTER_H__

#include <phool/PHObject.h>
#include "TrackerDefs.h"

#include <limits.h>
#include <cmath>
#include <iostream>
#include <set>

class TrackerCluster : public PHObject
{
 public:
  typedef std::set<TrackerDefs::keytype> HitSet;
  typedef std::set<TrackerDefs::keytype>::const_iterator ConstHitIter;
  typedef std::set<TrackerDefs::keytype>::iterator HitIter;

  virtual ~TrackerCluster() {}
  // PHObject virtual overloads

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrackerCluster base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual TrackerCluster* Clone() const { return NULL; }
  //
  // cluster id
  //
  virtual TrackerDefs::keytype get_id() const { return TrackerDefs::KEYMAX; }
  virtual void set_id(TrackerDefs::keytype id) {}
  //
  // cluster position
  //
  virtual float get_x() const { return NAN; }
  virtual void set_x(float x) {}
  virtual float get_y() const { return NAN; }
  virtual void set_y(float y) {}
  virtual float get_z() const { return NAN; }
  virtual void set_z(float z) {}
  virtual float get_position(int coor) const { return NAN; }
  virtual void set_position(int coor, float xi) {}
  virtual void set_global() {}
  virtual void set_local() {}
  virtual bool is_global() { return true; }
  //
  // cluster info
  //
  virtual float get_e() const { return NAN; }
  virtual void set_e(float e) {}
  virtual unsigned int get_adc() const { return UINT_MAX; }
  virtual void set_adc(unsigned int adc) {}
  virtual float get_size(unsigned int i, unsigned int j) const { return NAN; }
  virtual void set_size(unsigned int i, unsigned int j, float value) {}
  virtual float get_error(unsigned int i, unsigned int j) const { return NAN; }
  virtual void set_error(unsigned int i, unsigned int j, float value) {}
  //
  // clustered hit ids methods
  //
  virtual void clear_hits() {}
  virtual bool empty_hits() { return true; }
  virtual size_t size_hits() { return 0; }
  virtual void insert_hit(TrackerDefs::keytype hit_id) {}
  virtual size_t erase_hit(TrackerDefs::keytype hit_id) { return 0; }
  virtual ConstHitIter begin_hits() const { return HitSet().end(); }
  virtual ConstHitIter find_hit(TrackerDefs::keytype hitid) const { return HitSet().end(); }
  virtual ConstHitIter end_hits() const { return HitSet().end(); }
  virtual HitIter begin_hits() { return HitSet().end(); }
  virtual HitIter find_hit(TrackerDefs::keytype hitid) { return HitSet().end(); }
  virtual HitIter end_hits() { return HitSet().end(); }
  // convenience interface

  virtual float get_phi_size() const { return NAN; }
  virtual float get_z_size() const { return NAN; }
  virtual float get_phi_error() const { return NAN; }
  virtual float get_rphi_error() const { return NAN; }
  virtual float get_z_error() const { return NAN; }
 protected:
  TrackerCluster() {}
  ClassDef(TrackerCluster, 1);
};

#endif
