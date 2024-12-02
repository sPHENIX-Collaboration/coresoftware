// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCINFO_ZDCINFOV1_H
#define ZDCINFO_ZDCINFOV1_H

#include "Zdcinfo.h"

#include <limits>

class Zdcinfov1 : public Zdcinfo
{
 public:
  Zdcinfov1() = default;
  ~Zdcinfov1() override = default;

  void Reset() override { *this = Zdcinfov1(); }
  int isValid() const override;

  void set_zdc_energy(int arm, float zdc_e) override;
  float get_zdc_energy(int arm) const override;
  void set_radius(int arm, float _r) override;
  float get_radius(int arm) const override;

 private:
  float m_zdc_e[2] = {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  float m_radius[2] = {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
  ClassDefOverride(Zdcinfov1, 1);
};

#endif
