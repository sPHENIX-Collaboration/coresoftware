#ifndef G4BBC_BBCVERTEX_H
#define G4BBC_BBCVERTEX_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

class BbcVertex : public PHObject
{
 public:
  ~BbcVertex() override {}

  // PHObject virtual overloads

  void identify(std::ostream& os = std::cout) const override { os << "BbcVertex base class" << std::endl; }
  PHObject* CloneMe() const override { return nullptr; }
  int isValid() const  override{ return 0; }

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

 protected:
  BbcVertex() {}

 private:
  ClassDefOverride(BbcVertex, 1);
};

#endif
