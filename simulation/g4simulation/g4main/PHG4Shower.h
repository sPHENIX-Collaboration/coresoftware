// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SHOWER_H
#define G4MAIN_PHG4SHOWER_H

#include "PHG4HitDefs.h"

#include <phool/PHObject.h>

#include <cmath>  // for NAN def
#include <iostream>
#include <map>
#include <set>

class PHG4Shower : public PHObject
{
 public:
  typedef std::set<int> ParticleIdSet;
  typedef ParticleIdSet::iterator ParticleIdIter;
  typedef ParticleIdSet::const_iterator ParticleIdConstIter;

  typedef std::set<int> VertexIdSet;
  typedef VertexIdSet::iterator VertexIdIter;
  typedef VertexIdSet::const_iterator VertexIdConstIter;

  typedef std::map<int, std::set<PHG4HitDefs::keytype> > HitIdMap;
  typedef HitIdMap::iterator HitIdIter;
  typedef HitIdMap::const_iterator HitIdConstIter;

  ~PHG4Shower() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "PHG4Shower base class" << std::endl; }
  PHG4Shower* CloneMe() const override { return nullptr; }
  void Reset() override {}
  int isValid() const override { return 0; }

  // shower info

  virtual int get_id() const { return 0; }
  virtual void set_id(int /*id*/) {}

  virtual int get_parent_particle_id() const { return 0; }
  virtual void set_parent_particle_id(int /*parent_particle_id*/) {}

  virtual int get_parent_shower_id() const { return 0; }
  virtual void set_parent_shower_id(int /*parent_shower_id*/) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_position(unsigned int /*coor*/) const { return NAN; }
  virtual void set_position(unsigned int /*coor*/, float /*xi*/) {}

  virtual float get_covar(unsigned int /*i*/, unsigned int /*j*/) const { return NAN; }
  virtual void set_covar(unsigned int /*i*/, unsigned int /*j*/, float /*entry*/) {}

  virtual unsigned int get_nhits(int /*volume*/) const { return 0; }
  virtual void set_nhits(int /*volume*/, unsigned int /*nhits*/) {}

  virtual double get_edep() const { return NAN; }
  virtual float get_edep(int /*volume*/) const { return NAN; }
  virtual void set_edep(int /*volume*/, float /*edep*/) {}

  virtual double get_eion() const { return NAN; }
  virtual float get_eion(int /*volume*/) const { return NAN; }
  virtual void set_eion(int /*volume*/, float /*eion*/) {}

  virtual float get_light_yield(int /*volume*/) const { return NAN; }
  virtual void set_light_yield(int /*volume*/, float /*light_yield*/) {}

  virtual float get_eh_ratio(int /*volume*/) const { return NAN; }
  virtual void set_eh_ratio(int /*volume*/, float /*eh_ratio*/) {}

  virtual bool empty_g4particle_id() const { return true; }
  virtual size_t size_g4particle_id() const { return 0; }
  virtual void add_g4particle_id(int /*id*/) {}
  virtual ParticleIdIter begin_g4particle_id();
  virtual ParticleIdConstIter begin_g4particle_id() const;
  virtual ParticleIdIter end_g4particle_id();
  virtual ParticleIdConstIter end_g4particle_id() const;
  virtual size_t remove_g4particle_id(int /*id*/) { return 0; }
  virtual void clear_g4particle_id() {}
  virtual const ParticleIdSet& g4particle_ids() const = 0;

  virtual bool empty_g4vertex_id() const { return true; }
  virtual size_t size_g4vertex_id() const { return 0; }
  virtual void add_g4vertex_id(int /*id*/) {}
  virtual VertexIdIter begin_g4vertex_id();
  virtual VertexIdConstIter begin_g4vertex_id() const;
  virtual VertexIdIter end_g4vertex_id();
  virtual VertexIdConstIter end_g4vertex_id() const;
  virtual size_t remove_g4vertex_id(int /*id*/) { return 0; }
  virtual void clear_g4vertex_id() {}
  virtual const VertexIdSet& g4vertex_ids() const = 0;

  virtual bool empty_g4hit_id() const { return true; }
  virtual size_t size_g4hit_id() const { return 0; }
  virtual void add_g4hit_id(int /*volume*/, PHG4HitDefs::keytype /*id*/) {}
  virtual HitIdIter begin_g4hit_id();
  virtual HitIdConstIter begin_g4hit_id() const;
  virtual HitIdIter find_g4hit_id(int /*volume*/);
  virtual HitIdConstIter find_g4hit_id(int /*volume*/) const;
  virtual HitIdIter end_g4hit_id();
  virtual HitIdConstIter end_g4hit_id() const;
  virtual size_t remove_g4hit_id(int /*volume*/, PHG4HitDefs::keytype /*id*/) { return 0; }
  virtual size_t remove_g4hit_volume(int /*volume*/) { return 0; }
  virtual void clear_g4hit_id() {}
  virtual const HitIdMap& g4hit_ids() const = 0;

 protected:
  PHG4Shower() {}

 private:
  ClassDefOverride(PHG4Shower, 1);
};


/**
 * Equality operators for PHG4TruthInfoContainer. Note that the comparison is
 * performed only on the publicly accessible members.
 */
///@{
inline bool operator==(const PHG4Shower& lhs, const PHG4Shower& rhs)
{
  return lhs.get_id()                 == rhs.get_id()                 &&
         lhs.get_parent_particle_id() == rhs.get_parent_particle_id() &&
         lhs.get_parent_shower_id()   == rhs.get_parent_shower_id()   &&
         lhs.get_x()                  == rhs.get_x()                  &&
         lhs.get_y()                  == rhs.get_y()                  &&
         lhs.get_z()                  == rhs.get_z()                  &&
         lhs.g4particle_ids()         == rhs.g4particle_ids()         &&
         lhs.g4vertex_ids()           == rhs.g4vertex_ids()           &&
         lhs.g4hit_ids()              == rhs.g4hit_ids();
}

inline bool operator!=(const PHG4Shower& lhs, const PHG4Shower& rhs) { return !(lhs == rhs); }
///@}

#endif
