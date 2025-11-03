// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_TRUTHVERTEX_H
#define GLOBALVERTEX_TRUTHVERTEX_H

#include "Vertex.h"

#include <cmath>
#include <iostream>

class TruthVertex : public Vertex
{
 public:
  ~TruthVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "TruthVertex base class" << std::endl; }
  PHObject* CloneMe() const override { return nullptr; }
  int isValid() const override { return 0; }

  // vertex info

  virtual unsigned int get_id() const override { return std::numeric_limits<unsigned int>::max(); }
  virtual void set_id(unsigned int) override {}

  virtual float get_t() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_t(float) override {}

  virtual float get_t_err() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_t_err(float) override {}

  virtual float get_z() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_z(float) override {}

  virtual float get_z_err() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_z_err(float) override {}

  virtual float get_x() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_x(float) override {}

  virtual float get_x_err() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_x_err(float) {}

  virtual float get_y() const override { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_y(float) override {}

  virtual float get_y_err() const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_y_err(float) {}

 protected:
  TruthVertex() {}

 private:
  ClassDefOverride(TruthVertex, 1);
};

#endif