#ifndef G4BBC_BBCVERTEXV2_H
#define G4BBC_BBCVERTEXV2_H

#include "BbcVertex.h"

#include <iostream>

class BbcVertexv2 : public BbcVertex
{
 public:
  BbcVertexv2();
  ~BbcVertexv2() override;

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = BbcVertexv2(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new BbcVertexv2(*this); }

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

  void set_bbc_ns(int iarm, int bbc_npmt, float bbc_q, float bbc_t) override
  {
    _bbc_ns_npmt[iarm] = bbc_npmt;
    _bbc_ns_q[iarm] = bbc_q;
    _bbc_ns_t[iarm] = bbc_t;
  }

  int get_bbc_npmt(int iarm) const override { return _bbc_ns_npmt[iarm]; }
  float get_bbc_q(int iarm) const override { return _bbc_ns_q[iarm]; }
  float get_bbc_t(int iarm) const override { return _bbc_ns_t[iarm]; }

 private:
  unsigned int _id;  //< unique identifier within container
  float _t;          //< collision time
  float _t_err;      //< collision time uncertainty
  float _z;          //< collision position z
  float _z_err;      //< collision position z uncertainty
  int _bbc_ns_npmt[2];
  float _bbc_ns_q[2];
  float _bbc_ns_t[2];

  ClassDefOverride(BbcVertexv2, 1);
};

#endif
