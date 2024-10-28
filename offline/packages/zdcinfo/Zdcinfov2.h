// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCINFO_ZDCINFOV2_H
#define ZDCINFO_ZDCINFOV2_H

#include "Zdcinfo.h"

#include <limits>

class Zdcinfov2 : public Zdcinfo
{
 public:
  Zdcinfov2() = default;
  ~Zdcinfov2() override = default;

  void Reset() override { *this = Zdcinfov2(); }
int isValid() const override;

  void set_zdc_energy(int arm, float zdc_e) override;
  float get_zdc_energy(int arm) const override;
  void set_radius(int arm, float _r) override;
  float get_radius(int arm) const override;
  void set_zvertex(float _z) override; 
  float get_zvertex() const override; 

 private:
  float m_zdc_e[2] = {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  float m_radius[2] = {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  float m_zvertex = std::numeric_limits<float>::quiet_NaN();
  ClassDefOverride(Zdcinfov2, 1);
};

#endif
