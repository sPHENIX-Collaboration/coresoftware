#include "Zdcinfov2.h"

#include <cmath>

int Zdcinfov2::isValid() const
{
  if (std::isfinite(m_zdc_e[0]) && std::isfinite(m_zdc_e[1]))
  {
    return 1;
  }
  return 0;
}


void Zdcinfov2::set_zdc_energy(int arm, float zdc_e)
{
  m_zdc_e[arm] = zdc_e;
}

float Zdcinfov2::get_zdc_energy(int arm) const
{
  return m_zdc_e[arm];
}

void Zdcinfov2::set_radius(int arm, float _r)
{
  m_radius[arm] = _r;
}

float Zdcinfov2::get_radius(int arm) const
{
  return m_radius[arm];
}

void Zdcinfov2::set_zvertex(float _z)
{
  m_zvertex = _z;
}

float Zdcinfov2::get_zvertex() const
{
  return m_zvertex;
}


