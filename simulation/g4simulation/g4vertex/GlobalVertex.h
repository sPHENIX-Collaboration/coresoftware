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

  virtual ~GlobalVertex() {}

  // PHObject virtual overloads

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "GlobalVertex base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual GlobalVertex* Clone() const { return nullptr; }

  // vertex info

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int id) {}

  virtual float get_t() const { return NAN; }
  virtual void set_t(float t) {}

  virtual float get_t_err() const { return NAN; }
  virtual void set_t_err(float t_err) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float x) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float y) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float z) {}

  virtual float get_chisq() const { return NAN; }
  virtual void set_chisq(float chisq) {}

  virtual unsigned int get_ndof() const { return 0xFFFFFFFF; }
  virtual void set_ndof(unsigned int ndof) {}

  virtual float get_position(unsigned int coor) const { return NAN; }
  virtual void set_position(unsigned int coor, float xi) {}

  virtual float get_error(unsigned int i, unsigned int j) const { return NAN; }
  virtual void set_error(unsigned int i, unsigned int j, float value) {}

  //
  // associated vertex ids methods
  //
  virtual bool empty_vtxids() const { return true; }
  virtual size_t size_vtxids() const { return 0; }
  virtual size_t count_vtxids(VTXTYPE type) const { return 0; }

  virtual void clear_vtxids() {}
  virtual void insert_vtxids(VTXTYPE type, unsigned int vtxid) {}
  virtual size_t erase_vtxids(VTXTYPE type) { return 0; }
  virtual void erase_vtxids(VtxIter iter) {}
  virtual void erase_vtxids(VtxIter first, VtxIter last) {}

  virtual ConstVtxIter begin_vtxids() const { return std::map<VTXTYPE, unsigned int>().end(); }
  virtual ConstVtxIter find_vtxids(VTXTYPE type) const { return std::map<VTXTYPE, unsigned int>().end(); }
  virtual ConstVtxIter end_vtxids() const { return std::map<VTXTYPE, unsigned int>().end(); }

  virtual VtxIter begin_vtxids() { return std::map<VTXTYPE, unsigned int>().end(); }
  virtual VtxIter find_vtxids(VTXTYPE type) { return std::map<VTXTYPE, unsigned int>().end(); }
  virtual VtxIter end_vtxids() { return std::map<VTXTYPE, unsigned int>().end(); }

 protected:
  GlobalVertex() {}

 private:
  ClassDef(GlobalVertex, 1);
};

#endif
