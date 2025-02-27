// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_GLOBALVERTEX_H
#define GLOBALVERTEX_GLOBALVERTEX_H

#include "Vertex.h"

#include <phool/PHObject.h>

#include <iostream>
#include <limits>
#include <map>

class GlobalVertex : public PHObject
{
 public:
  // the order matters (best vertex -> highest number), so leave some space in case we want to wedge other vertices in here
  enum VTXTYPE
  {
    UNDEFINED = 0,
    TRUTH = 100,
    SMEARED = 200,
    MBD = 300,
    SVTX = 400,
    SVTX_MBD = 500
  };

  typedef std::vector<const Vertex*> VertexVector;
  typedef std::map<GlobalVertex::VTXTYPE, VertexVector>::const_iterator ConstVertexIter;
  typedef std::map<GlobalVertex::VTXTYPE, VertexVector>::iterator VertexIter;

  // Deprecated as of GlobalVertexv2
  typedef std::map<GlobalVertex::VTXTYPE, unsigned int>::const_iterator ConstVtxIter;
  typedef std::map<GlobalVertex::VTXTYPE, unsigned int>::iterator VtxIter;

  ~GlobalVertex() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override
  {
    os << "GlobalVertex base class" << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  // vertex info

  virtual unsigned int get_id() const { return std::numeric_limits<unsigned int>::max(); }
  virtual void set_id(unsigned int) { return; }

  virtual float get_t() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_t(float) { return; }

  virtual float get_t_err() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_t_err(float) { return; }

  virtual float get_x() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_x(float) { return; }

  virtual float get_y() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_y(float) { return; }

  virtual float get_z() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_z(float) { return; }

  virtual float get_chisq() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_chisq(float) { return; }

  virtual unsigned int get_ndof() const { return std::numeric_limits<unsigned int>::max(); }
  virtual void set_ndof(unsigned int) { return; }

  virtual float get_position(unsigned int /*coor*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_position(unsigned int /*coor*/, float /*xi*/) { return; }

  virtual float get_error(unsigned int /*i*/, unsigned int /*j*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_error(unsigned int /*i*/, unsigned int /*j*/, float /*value*/) { return; }

  virtual unsigned int get_beam_crossing() const { return std::numeric_limits<unsigned int>::max(); }
  virtual void set_beam_crossing(unsigned int) { return; }

  virtual bool empty_vtxs() const { return true; }
  virtual size_t size_vtxs() const { return 0; }
  virtual size_t count_vtxs(VTXTYPE) const { return 0; }
  virtual void clear_vtxs() { return; }
  virtual void insert_vtx(VTXTYPE, const Vertex*) { return; }
  virtual void clone_insert_vtx(VTXTYPE, const Vertex*) { return; }
  virtual size_t erase_vtxs(VTXTYPE) { return 0; }
  virtual void erase_vtxs(VertexIter) { return; }

  virtual ConstVertexIter begin_vertexes() const;
  virtual ConstVertexIter find_vertexes(VTXTYPE type) const;
  virtual ConstVertexIter end_vertexes() const;

  virtual VertexIter begin_vertexes();
  virtual VertexIter find_vertexes(VTXTYPE type);
  virtual VertexIter end_vertexes();

  //
  // associated vertex ids methods
  // vtx id container accessors are deprecated.
  // Use actual vertex container accessors instead
  //
  virtual bool empty_vtxids() const { return true; }
  virtual size_t size_vtxids() const { return 0; }
  virtual size_t count_vtxids(VTXTYPE /*type*/) const { return 0; }

  virtual void clear_vtxids() { return; }
  virtual void insert_vtxids(VTXTYPE /*type*/, unsigned int /*vtxid*/) { return; }
  virtual size_t erase_vtxids(VTXTYPE /*type*/) { return 0; }
  virtual void erase_vtxids(VtxIter /*iter*/) { return; }
  virtual void erase_vtxids(VtxIter /*first*/, VtxIter /*last*/) { return; }

  virtual ConstVtxIter begin_vtxids() const;
  virtual ConstVtxIter find_vtxids(VTXTYPE type) const;
  virtual ConstVtxIter end_vtxids() const;

  virtual VtxIter begin_vtxids();
  virtual VtxIter find_vtxids(VTXTYPE type);
  virtual VtxIter end_vtxids();

 protected:
  GlobalVertex() = default;

 private:
  ClassDefOverride(GlobalVertex, 1);
};

#endif
