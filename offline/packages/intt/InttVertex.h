#ifndef INTT_INTTVERTEX_H
#define INTT_INTTVERTEX_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

class InttVertex : public PHObject
{
 public:
  ~InttVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "InttVertex base class" << std::endl; }
  PHObject* CloneMe() const override { return nullptr; }
  int isValid() const override { return 0; }

  // vertex info

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_z_err() const { return NAN; }
  virtual void set_z_err(float) {}

 protected:
  InttVertex() {}

 private:
  ClassDefOverride(InttVertex, 1);
};

#endif
