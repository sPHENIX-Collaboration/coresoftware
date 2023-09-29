#ifndef INTT_INTTVERTEXV1_H
#define INTT_INTTVERTEXV1_H

#include "InttVertex.h"

#include <iostream>

class PHObject;

class InttVertexv1 : public InttVertex
{
 public:
  InttVertexv1() = default;
  ~InttVertexv1() override = default;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = InttVertexv1(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new InttVertexv1(*this); }

  // vertex info

  unsigned int get_id() const override { return _id; }
  void set_id(unsigned int id) override { _id = id; }

  float get_z() const override { return _z; }
  void set_z(float z) override { _z = z; }

  float get_z_err() const override { return _z_err; }
  void set_z_err(float z_err) override { _z_err = z_err; }

 private:
  unsigned int _id = std::numeric_limits<unsigned int>::max();  //< unique identifier within container
  float _z = std::numeric_limits<float>::signaling_NaN();          //< collision position z
  float _z_err = std::numeric_limits<float>::signaling_NaN();      //< collision position z uncertainty

  ClassDefOverride(InttVertexv1, 1);
};

#endif
