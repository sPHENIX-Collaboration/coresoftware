#include "EpdGeomV1.h"

#include <calobase/TowerInfoDefs.h>

#include <phool/PHObject.h>

#include <map>
#include <utility>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>


float EpdGeomV1::get_r(unsigned int key) const
{
  return tile_r[TowerInfoDefs::get_epd_rbin(key)];
}

float EpdGeomV1::get_z(unsigned int key) const
{
  return tile_z[TowerInfoDefs::get_epd_arm(key)];
}

float EpdGeomV1::get_phi(unsigned int key) const
{

  if(TowerInfoDefs::get_epd_rbin(key) == 0)
  {
    return tile_phi0[TowerInfoDefs::get_epd_phibin(key)];
  }
  else
  {
    return tile_phi[TowerInfoDefs::get_epd_phibin(key)];
  }

}
 
void EpdGeomV1::set_z(unsigned int key, float z)
{
  tile_z[TowerInfoDefs::get_epd_arm(key)] = z;
}

void EpdGeomV1::set_r(unsigned int key, float r)
{
  tile_r[TowerInfoDefs::get_epd_rbin(key)] = r;
}

void EpdGeomV1::set_phi(unsigned int key, float f)
{
  tile_phi[TowerInfoDefs::get_epd_phibin(key)] = f;
}

void EpdGeomV1::set_phi0(unsigned int key, float f0)
{
  tile_phi0[TowerInfoDefs::get_epd_phibin(key)] = f0;
}

