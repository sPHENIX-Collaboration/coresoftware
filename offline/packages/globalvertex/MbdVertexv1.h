#ifndef G4MBD_MBDVERTEXV1_H
#define G4MBD_MBDVERTEXV1_H

#include "MbdVertex.h"

#include <iostream>

class MbdVertexv1 : public MbdVertex
{
 public:
  MbdVertexv1();
  ~MbdVertexv1() override;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = MbdVertexv1(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new MbdVertexv1(*this); }

  // vertex info

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  float get_t() const override { return _t; }
  void set_t(float t) override { _t = t; }

  float get_t_err() const override { return _t_err; }
  void set_t_err(float t_err) override { _t_err = t_err; }

  float get_z() const override { return _z; }
  void set_z(float z) override { _z = z; }

  float get_z_err() const override { return _z_err; }
  void set_z_err(float z_err) override { _z_err = z_err; }

 private:
  unsigned int _id;  //< unique identifier within container
  float _t;          //< collision time
  float _t_err;      //< collision time uncertainty
  float _z;          //< collision position z
  float _z_err;      //< collision position z uncertainty

  ClassDefOverride(MbdVertexv1, 1);
};

#endif
