#ifndef G4MBD_MBDVERTEXV2_H
#define G4MBD_MBDVERTEXV2_H

#include "MbdVertex.h"

#include <iostream>

class MbdVertexv2 : public MbdVertex
{
 public:
  MbdVertexv2();
  ~MbdVertexv2() override;

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

  unsigned int get_beam_crossing() const override { return _bco; }
  void set_beam_crossing(unsigned int bco) override { _bco = bco; }

 private:
  unsigned int _id;   //< unique identifier within container
  unsigned int _bco;  //< global bco
  float _t;           //< collision time
  float _t_err;       //< collision time uncertainty
  float _z;           //< collision position z
  float _z_err;       //< collision position z uncertainty

  ClassDefOverride(MbdVertexv2, 1);
};

#endif
