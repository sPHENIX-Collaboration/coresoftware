#ifndef __SVTXCLUSTER_H__
#define __SVTXCLUSTER_H__

#include <phool/PHObject.h>

#include <cmath>
#include <climits>
#include <iostream>
#include <set>

class SvtxCluster : public PHObject
{
 public:
  typedef std::set<unsigned int> HitSet;
  typedef std::set<unsigned int>::const_iterator ConstHitIter;
  typedef std::set<unsigned int>::iterator HitIter;

  virtual ~SvtxCluster() {}

  // PHObject virtual overloads

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "SvtxCluster base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual PHObject* CloneMe() const { return nullptr; }

  // cluster info

  virtual unsigned int get_id() const { return UINT_MAX; }
  virtual void set_id(unsigned int id) {}

  virtual unsigned int get_layer() const { return UINT_MAX; }
  virtual void set_layer(unsigned int layer) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float x) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float y) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float z) {}

  virtual float get_position(int coor) const { return NAN; }
  virtual void set_position(int coor, float xi) {}

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
  virtual void insert_hit(unsigned int hit_id) {}
  virtual size_t erase_hit(unsigned int hit_id) { return 0; }
  virtual ConstHitIter begin_hits() const { return HitSet().end(); }
  virtual ConstHitIter find_hit(unsigned int hitid) const { return HitSet().end(); }
  virtual ConstHitIter end_hits() const { return HitSet().end(); }
  virtual HitIter begin_hits() { return HitSet().end(); }
  virtual HitIter find_hit(unsigned int hitid) { return HitSet().end(); }
  virtual HitIter end_hits() { return HitSet().end(); }

  // convenience interface

  virtual float get_phi_size() const { return NAN; }
  virtual float get_z_size() const { return NAN; }

  virtual float get_phi_error() const { return NAN; }
  virtual float get_rphi_error() const { return NAN; }
  virtual float get_z_error() const { return NAN; }

 protected:
  SvtxCluster() {}

  ClassDef(SvtxCluster, 1);
};

#endif
