#ifndef __TrkrCluster_H__
#define __TrkrCluster_H__

#include <phool/PHObject.h>
#include "TrkrDefUtil.h"

#include <limits.h>
#include <cmath>
#include <iostream>
#include <set>

class TrkrCluster : public PHObject
{
 public:
  typedef std::set<TrkrDefs::hitsetkey> HitSet;
  typedef std::set<TrkrDefs::hitsetkey>::const_iterator ConstHitIter;
  typedef std::set<TrkrDefs::hitsetkey>::iterator HitIter;

  virtual ~TrkrCluster() {}
  // PHObject virtual overloads

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrCluster base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual TrkrCluster* Clone() const { return NULL; }
  //
  // cluster id
  //
  virtual TrkrDefs::cluskey get_id() const { return TrkrDefs::CLUSKEYMAX; }
  virtual void set_id(TrkrDefs::cluskey id) {}
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
  virtual void insert_hit(TrkrDefs::hitsetkey hit_id) {}
  virtual size_t erase_hit(TrkrDefs::hitsetkey hit_id) { return 0; }
  virtual ConstHitIter begin_hits() const { return HitSet().end(); }
  virtual ConstHitIter find_hit(TrkrDefs::hitsetkey hitid) const { return HitSet().end(); }
  virtual ConstHitIter end_hits() const { return HitSet().end(); }
  virtual HitIter begin_hits() { return HitSet().end(); }
  virtual HitIter find_hit(TrkrDefs::hitsetkey hitid) { return HitSet().end(); }
  virtual HitIter end_hits() { return HitSet().end(); }
  //
  // convenience interface
  //
  virtual float get_phi_size() const { return NAN; }
  virtual float get_z_size() const { return NAN; }
  virtual float get_phi_error() const { return NAN; }
  virtual float get_rphi_error() const { return NAN; }
  virtual float get_z_error() const { return NAN; }
 protected:
  TrkrCluster() {}
  ClassDef(TrkrCluster, 1);
};

#endif
