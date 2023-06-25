#include "RawTowerGeomv2.h"

#include <cmath>
#include <iostream>

RawTowerGeomv2::RawTowerGeomv2(RawTowerDefs::keytype id)
  : _towerid(id)
{
}

double RawTowerGeomv2::get_center_radius() const
{
  return std::sqrt(_center_x * _center_x +
                   _center_y * _center_y);
}

double RawTowerGeomv2::get_eta() const
{
  double eta;
  double radius;
  double theta;
  radius = std::sqrt(_center_x * _center_x + _center_y * _center_y);
  theta = std::atan2(radius, _center_z);
  eta = -std::log(std::tan(theta / 2.));

  return eta;
}

double RawTowerGeomv2::get_phi() const
{
  return std::atan2(_center_y, _center_x);
}

void RawTowerGeomv2::identify(std::ostream& os) const
{
  os << "RawTowerGeomv2:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z()
     << "\n           dx: " << get_size_x() << " dy: " << get_size_y() << " dz: " << get_size_z() << std::endl;
}
