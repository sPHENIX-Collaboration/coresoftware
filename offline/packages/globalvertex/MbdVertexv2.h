// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_MBDVERTEXV2_H
#define GLOBALVERTEX_MBDVERTEXV2_H

#include "MbdVertex.h"

#include <iostream>
#include <limits>

class MbdVertexv2 : public MbdVertex
{
 public:
  MbdVertexv2() = default;
  ~MbdVertexv2() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = MbdVertexv2(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new MbdVertexv2(*this); }

  // vertex info

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  float get_t() const override { return _t; }
  void set_t(float t) override { _t = t; }

  float get_t_err() const override { return _t_err; }
  void set_t_err(float t_err) override { _t_err = t_err; }

  // Return 0 for now, can implement beam spot
  float get_x() const override { return 0; }
  float get_y() const override { return 0; }

  float get_z() const override { return _z; }
  void set_z(float z) override { _z = z; }

  float get_z_err() const override { return _z_err; }
  void set_z_err(float z_err) override { _z_err = z_err; }

  float get_position(unsigned int coor) const override;

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

 private:
  static constexpr short int short_int_max = std::numeric_limits<short int>::max();  // 32767

  static short int rollover_short(short int bco)
  {
    if (bco == short_int_max) return short_int_max;
    if (bco >= 0) return bco;

    const int bco_ro = static_cast<int>(short_int_max) + static_cast<int>(bco);  // bco negative
    return static_cast<short int>(bco_ro);
  }

  static short int rollover_from_unsignedint(unsigned int bco)
  {
    // if unsigned int max, return short int max
    if (bco == std::numeric_limits<unsigned int>::max())
    {
      return short_int_max;
    }

    // common case: [0, 32767]
    if (bco <= static_cast<unsigned int>(short_int_max))
    {
      return static_cast<short int>(bco);
    }

    const short int bco_ro = static_cast<short int>(static_cast<unsigned short>(bco));
    if (bco_ro >= 0) return bco_ro;

    return rollover_short(bco_ro);
  }

  unsigned int _id{std::numeric_limits<unsigned int>::max()};   //< unique identifier within container
  unsigned int _bco{std::numeric_limits<unsigned int>::max()};  //< global bco (legacy storage)
  float _t{std::numeric_limits<float>::quiet_NaN()};            //< collision time
  float _t_err{std::numeric_limits<float>::quiet_NaN()};        //< collision time uncertainty
  float _z{std::numeric_limits<float>::quiet_NaN()};            //< collision position z
  float _z_err{std::numeric_limits<float>::quiet_NaN()};        //< collision position z uncertainty

  ClassDefOverride(MbdVertexv2, 1);
};

#endif
