// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCINFO_ZDCINFO_H
#define ZDCINFO_ZDCINFO_H

#include <phool/PHObject.h>

#include <limits>

class Zdcinfo : public PHObject
{
 public:
  ~Zdcinfo() override {}

  virtual void set_zdc_energy(int /*arm*/, float /*zdc_e*/) { return; }
  virtual float get_zdc_energy(const int /*arm*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_radius(int /*arm*/, float /*_r*/) { return; }
  virtual float get_radius(const int /*arm*/) const { return std::numeric_limits<float>::quiet_NaN(); }
  virtual void set_zvertex(float /*_z*/) { return; }
  virtual float get_zvertex() const { return std::numeric_limits<float>::quiet_NaN(); }
 
 protected:
  Zdcinfo() {}

 private:
  ClassDefOverride(Zdcinfo, 1);
};

#endif
