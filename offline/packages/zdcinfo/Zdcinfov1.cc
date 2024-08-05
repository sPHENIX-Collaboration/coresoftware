#include "Zdcinfov1.h"
#include <cmath>

void Zdcinfov1::set_zdc_energy(int arm, float zdc_e)
{
  m_zdc_e[arm] = zdc_e;
}

float Zdcinfov1::get_zdc_energy(int arm) const
{
  return m_zdc_e[arm];
}

void Zdcinfov1::set_radius(int arm, float _r)
{
  m_radius[arm] = _r;
}

float Zdcinfov1::get_radius(int arm) const
{
  return m_radius[arm];
}
