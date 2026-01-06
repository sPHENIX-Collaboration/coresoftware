#include "EpdGeomV2.h"

#include <calobase/TowerInfoDefs.h>

#include <phool/PHObject.h>

#include <array>

float EpdGeomV2::get_r(unsigned int key) const
{
  return tile_r.at(TowerInfoDefs::get_epd_rbin(key));
}

float EpdGeomV2::get_z(unsigned int key) const
{
  return tile_z.at(TowerInfoDefs::get_epd_arm(key));
}

float EpdGeomV2::get_phi(unsigned int key) const
{
  if (TowerInfoDefs::get_epd_rbin(key) == 0)
  {
    return tile_phi0.at(TowerInfoDefs::get_epd_phibin(key));
  }

  return tile_phi.at(TowerInfoDefs::get_epd_phibin(key));
}

void EpdGeomV2::set_z(unsigned int key, float z)
{
  tile_z.at(TowerInfoDefs::get_epd_arm(key)) = z;
}

void EpdGeomV2::set_r(unsigned int key, float r)
{
  tile_r.at(TowerInfoDefs::get_epd_rbin(key)) = r;
}

void EpdGeomV2::set_phi(unsigned int key, float f)
{
  tile_phi.at(TowerInfoDefs::get_epd_phibin(key)) = f;
}

void EpdGeomV2::set_phi0(unsigned int key, float f0)
{
  tile_phi0.at(TowerInfoDefs::get_epd_phibin(key)) = f0;
}
