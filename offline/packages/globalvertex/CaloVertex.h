// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_CALOVERTEX_H
#define GLOBALVERTEX_CALOVERTEX_H

#include "Vertex.h"

#include <cmath>
#include <iostream>

class CaloVertex : public Vertex
{
 public:
  ~CaloVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "CaloVertex base class" << std::endl; }
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

 protected:
  CaloVertex() {}

 private:
  ClassDefOverride(CaloVertex, 1);
};

#endif
