// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_GLOBALVERTEXV3_H
#define GLOBALVERTEX_GLOBALVERTEXV3_H

#include "GlobalVertex.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <map>

class PHObject;

class GlobalVertexv3 : public GlobalVertex
{
 public:
  GlobalVertexv3() = default;
  GlobalVertexv3(const unsigned int id);
  ~GlobalVertexv3() override;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new GlobalVertexv3(*this); }

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  short int get_beam_crossing() const override { return _bco; }
  void set_beam_crossing(short int bco) override { _bco = bco; }

  float get_t() const override;
  float get_t_err() const override;
  float get_x() const override;
  float get_y() const override;
  float get_z() const override;
  float get_chisq() const override;
  unsigned int get_ndof() const override;
  float get_position(unsigned int coor) const override;
  float get_error(unsigned int i, unsigned int j) const override;

  //
  // associated vertex methods
  //
  bool empty_vtxs() const override { return _vtxs.empty(); }
  size_t size_vtxs() const override { return _vtxs.size(); }
  size_t count_vtxs(GlobalVertex::VTXTYPE type) const override;

  void clear_vtxs() override { _vtxs.clear(); }
  void insert_vtx(GlobalVertex::VTXTYPE type, const Vertex* vertex) override;
  void clone_insert_vtx(GlobalVertex::VTXTYPE type, const Vertex* vertex) override;
  size_t erase_vtxs(GlobalVertex::VTXTYPE type) override { return _vtxs.erase(type); }
  void erase_vtxs(GlobalVertex::VertexIter iter) override { _vtxs.erase(iter); }

  GlobalVertex::ConstVertexIter begin_vertexes() const override { return _vtxs.begin(); }
  GlobalVertex::ConstVertexIter find_vertexes(GlobalVertex::VTXTYPE type) const override { return _vtxs.find(type); }
  GlobalVertex::ConstVertexIter end_vertexes() const override { return _vtxs.end(); }

  GlobalVertex::VertexIter begin_vertexes() override { return _vtxs.begin(); }
  GlobalVertex::VertexIter find_vertexes(GlobalVertex::VTXTYPE type) override { return _vtxs.find(type); }
  GlobalVertex::VertexIter end_vertexes() override { return _vtxs.end(); }

 private:
  unsigned int _id{std::numeric_limits<unsigned int>::max()};
  short int _bco{std::numeric_limits<short int>::max()};          //< global bco (signed short)
  std::map<GlobalVertex::VTXTYPE, VertexVector> _vtxs;            //< list of vtxs

  ClassDefOverride(GlobalVertexv3, 3);
};

#endif
