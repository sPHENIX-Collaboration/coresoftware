// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_GLOBALVERTEXV2_H
#define GLOBALVERTEX_GLOBALVERTEXV2_H

#include "GlobalVertex.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <limits>
#include <map>

class PHObject;

class GlobalVertexv2 : public GlobalVertex
{
 public:
  GlobalVertexv2() = default;
  GlobalVertexv2(const unsigned int id);
  ~GlobalVertexv2() override;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override;
  PHObject* CloneMe() const override { return new GlobalVertexv2(*this); }

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  short int get_beam_crossing() const override
  {
    return rollover_from_unsignedint(_bco);
  }

  void set_beam_crossing(short int bco) override
  {
    if (bco == short_int_max)
    {
      _bco = std::numeric_limits<unsigned int>::max();
      return;
    }

    const short int bco_ro = rollover_short(bco);
    _bco = static_cast<unsigned int>(bco_ro);
  }

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
  static constexpr short int short_int_max = std::numeric_limits<short int>::max();

  static short int rollover_short(short int bco)
  {
    if (bco == short_int_max) return short_int_max;
    if (bco >= 0) return bco;
    return static_cast<short int>(static_cast<int>(short_int_max) + static_cast<int>(bco));
  }

  static short int rollover_from_unsignedint(unsigned int bco)
  {
    if (bco == std::numeric_limits<unsigned int>::max()) return short_int_max;
    if (bco <= static_cast<unsigned int>(short_int_max)) return static_cast<short int>(bco);

    const short int bco_ro = static_cast<short int>(static_cast<unsigned short>(bco));
    if (bco_ro >= 0) return bco_ro;
    return rollover_short(bco_ro);
  }

  unsigned int _id{std::numeric_limits<unsigned int>::max()};
  unsigned int _bco{std::numeric_limits<unsigned int>::max()};  //< global bco
  std::map<GlobalVertex::VTXTYPE, VertexVector> _vtxs;          //< list of vtxs

  ClassDefOverride(GlobalVertexv2, 2);
};

#endif
