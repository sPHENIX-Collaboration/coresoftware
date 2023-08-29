#ifndef INTT_INTTVERTEXV1_H
#define INTT_INTTVERTEXV1_H

#include "InttVertex.h"

#include <iostream>

class InttVertexv1 : public InttVertex
{
 public:
  InttVertexv1();
  ~InttVertexv1() override;

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
  unsigned int _id;  //< unique identifier within container
  float _z;          //< collision position z
  float _z_err;      //< collision position z uncertainty

  ClassDefOverride(InttVertexv1, 1);
};

#endif
