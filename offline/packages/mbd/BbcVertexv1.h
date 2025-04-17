#ifndef G4BBC_BBCVERTEXV1_H
#define G4BBC_BBCVERTEXV1_H

#include "BbcVertex.h"

#include <iostream>
#include <limits>

class BbcVertexv1 : public BbcVertex
{
 public:
  BbcVertexv1() = default;;
  ~BbcVertexv1() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = BbcVertexv1(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new BbcVertexv1(*this); }

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
  unsigned int _id {std::numeric_limits<unsigned int>::max()};  //< unique identifier within container
  float _t {std::numeric_limits<float>::quiet_NaN()};          //< collision time
  float _t_err {std::numeric_limits<float>::quiet_NaN()};      //< collision time uncertainty
  float _z {std::numeric_limits<float>::quiet_NaN()};          //< collision position z
  float _z_err {std::numeric_limits<float>::quiet_NaN()};      //< collision position z uncertainty

  ClassDefOverride(BbcVertexv1, 1);
};

#endif
