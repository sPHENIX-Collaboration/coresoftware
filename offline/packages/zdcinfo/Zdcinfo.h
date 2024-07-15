// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCINFO_H
#define ZDCINFO_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>

class Zdcinfo : public PHObject
{
 public:
  ~Zdcinfo() override {}

  virtual void set_zdc_energy(int /*arm*/, float /*zdc_e*/) { return; }
  virtual float get_zdc_energy(const int /*arm*/) const { return NAN; }
  virtual void set_radius(int /*arm*/, float /*_r*/) { return; }
  virtual float get_radius(const int /*arm*/) const { return NAN; }

 protected:
  Zdcinfo() {}

 private:
  ClassDefOverride(Zdcinfo, 1);
};

#endif
