#ifndef GLOBALVERTEX_TRUTHVERTEXV1_H
#define GLOBALVERTEX_TRUTHVERTEXV1_H

#include "TruthVertex.h"

#include <iostream>
#include <limits>

class TruthVertex_v1 : public TruthVertex
{
 public:
  TruthVertex_v1() = default;
  ~TruthVertex_v1() override = default;

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = TruthVertex_v1(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new TruthVertex_v1(*this); }

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
  unsigned int _id{std::numeric_limits<unsigned int>::max()};
  float _t{std::numeric_limits<float>::quiet_NaN()};
  float _t_err{std::numeric_limits<float>::quiet_NaN()};
  float _z{std::numeric_limits<float>::quiet_NaN()};
  float _z_err{std::numeric_limits<float>::quiet_NaN()};

  ClassDefOverride(TruthVertex_v1, 1);
};

#endif
