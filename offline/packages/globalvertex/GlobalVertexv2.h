// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4VERTEX_GLOBALVERTEXV2_H
#define G4VERTEX_GLOBALVERTEXV2_H

#include "GlobalVertex.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <map>
#include <utility>  // for pair, make_pair

class PHObject;

class GlobalVertexv2 : public GlobalVertex
{
 public:
  GlobalVertexv2();
  GlobalVertexv2(const unsigned int id);
  ~GlobalVertexv2() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = GlobalVertexv2(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new GlobalVertexv2(*this); }

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  unsigned int get_beam_crossing() const override { return _bco; }
  void set_beam_crossing(unsigned int bco) override { _bco = bco; }

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
  size_t erase_vtxs(GlobalVertex::VTXTYPE type) override { return _vtxs.erase(type); }
  void erase_vtxs(GlobalVertex::VertexIter iter) override { _vtxs.erase(iter); }

  GlobalVertex::ConstVertexIter begin_vertexes() const override { return _vtxs.begin(); }
  GlobalVertex::ConstVertexIter find_vertexes(GlobalVertex::VTXTYPE type) const override { return _vtxs.find(type); }
  GlobalVertex::ConstVertexIter end_vertexes() const override { return _vtxs.end(); }

  GlobalVertex::VertexIter begin_vertexes() override { return _vtxs.begin(); }
  GlobalVertex::VertexIter find_vertexes(GlobalVertex::VTXTYPE type) override { return _vtxs.find(type); }
  GlobalVertex::VertexIter end_vertexes() override { return _vtxs.end(); }

 private:
  unsigned int _id;
  unsigned int _bco;  //< global bco
  std::map<GlobalVertex::VTXTYPE, VertexVector> _vtxs;  //< list of vtxs

  ClassDefOverride(GlobalVertexv2, 2);
};

#endif
