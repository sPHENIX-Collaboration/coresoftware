// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4VERTEX_GLOBALVERTEX_H
#define G4VERTEX_GLOBALVERTEX_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <map>

class GlobalVertex : public PHObject
{
 public:
  enum VTXTYPE
  {
    NONE = 0,
    BBC = 1,
    SVTX = 2
  };

  typedef std::map<GlobalVertex::VTXTYPE, unsigned int>::const_iterator ConstVtxIter;
  typedef std::map<GlobalVertex::VTXTYPE, unsigned int>::iterator VtxIter;

  ~GlobalVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override
  {
    os << "GlobalVertex base class" << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  // vertex info

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) {}

  virtual float get_t() const { return NAN; }
  virtual void set_t(float) {}

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

  virtual unsigned int get_ndof() const { return 0xFFFFFFFF; }
  virtual void set_ndof(unsigned int) {}

  virtual float get_position(unsigned int /*coor*/) const { return NAN; }
  virtual void set_position(unsigned int /*coor*/, float /*xi*/) {}

  virtual float get_error(unsigned int /*i*/, unsigned int /*j*/) const { return NAN; }
  virtual void set_error(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) {}

  //
  // associated vertex ids methods
  //
  virtual bool empty_vtxids() const { return true; }
  virtual size_t size_vtxids() const { return 0; }
  virtual size_t count_vtxids(VTXTYPE /*type*/) const { return 0; }

  virtual void clear_vtxids() {}
  virtual void insert_vtxids(VTXTYPE /*type*/, unsigned int /*vtxid*/) {}
  virtual size_t erase_vtxids(VTXTYPE /*type*/) { return 0; }
  virtual void erase_vtxids(VtxIter /*iter*/) {}
  virtual void erase_vtxids(VtxIter /*first*/, VtxIter /*last*/) {}

  virtual ConstVtxIter begin_vtxids() const;
  virtual ConstVtxIter find_vtxids(VTXTYPE type) const;
  virtual ConstVtxIter end_vtxids() const;

  virtual VtxIter begin_vtxids();
  virtual VtxIter find_vtxids(VTXTYPE type);
  virtual VtxIter end_vtxids();

 protected:
  GlobalVertex() {}

 private:
  ClassDefOverride(GlobalVertex, 1);
};

#endif
