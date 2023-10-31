#ifndef G4MBD_MBDVERTEX_H
#define G4MBD_MBDVERTEX_H

#include "Vertex.h"

#include <cmath>
#include <iostream>

class MbdVertex : public Vertex
{
 public:
  ~MbdVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "MbdVertex base class" << std::endl; }
  PHObject* CloneMe() const override { return nullptr; }
  int isValid() const override { return 0; }

  // vertex info

  virtual unsigned int get_id() const override { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) override {}

  virtual float get_t() const override { return NAN; }
  virtual void set_t(float) override {}

  virtual float get_t_err() const override { return NAN; }
  virtual void set_t_err(float) override {}

  virtual float get_z() const override { return NAN; }
  virtual void set_z(float) override {}

  virtual float get_z_err() const override { return NAN; }
  virtual void set_z_err(float) override {}

  virtual unsigned int get_beam_crossing() const override { return std::numeric_limits<unsigned int>::max(); }
  virtual void set_beam_crossing(unsigned int) override {}

  virtual void set_bbc_ns(int, int, float, float) override {}
  virtual int get_bbc_npmt(int) const override { return std::numeric_limits<int>::max(); }
  virtual float get_bbc_q(int) const override { return NAN; }
  virtual float get_bbc_t(int) const override { return NAN; }

 protected:
  MbdVertex() {}

 private:
  ClassDefOverride(MbdVertex, 1);
};

#endif
