#ifndef G4MBD_MBDVERTEX_H
#define G4MBD_MBDVERTEX_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

class MbdVertex : public PHObject
{
 public:
  ~MbdVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "MbdVertex base class" << std::endl; }
  PHObject* CloneMe() const override { return nullptr; }
  int isValid() const override { return 0; }

  // vertex info

  virtual unsigned int get_id() const { return 0xFFFFFFFF; }
  virtual void set_id(unsigned int) {}

  virtual float get_t() const { return NAN; }
  virtual void set_t(float) {}

  virtual float get_t_err() const { return NAN; }
  virtual void set_t_err(float) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_z_err() const { return NAN; }
  virtual void set_z_err(float) {}

  virtual void set_bbc_ns(int, int, float, float) {}
  virtual int get_bbc_npmt(int) const { return std::numeric_limits<int>::max(); }
  virtual float get_bbc_q(int) const { return NAN; }
  virtual float get_bbc_t(int) const { return NAN; }

 protected:
  MbdVertex() {}

 private:
  ClassDefOverride(MbdVertex, 1);
};

#endif
