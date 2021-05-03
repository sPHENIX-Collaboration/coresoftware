#include "RawTowerGeomv1.h"

#include <cmath>
#include <iostream>

RawTowerGeomv1::RawTowerGeomv1(RawTowerDefs::keytype id)
  : _towerid(id)
{
}

double RawTowerGeomv1::get_center_radius() const
{
  return std::sqrt(_center_x * _center_x +
                   _center_y * _center_y);
}

double RawTowerGeomv1::get_theta() const
{
  double radius = sqrt(_center_x * _center_x + _center_y * _center_y);
  double theta = atan2(radius, _center_z);
  return theta;
}

double RawTowerGeomv1::get_eta() const
{
  double theta = get_theta();
  double eta = -log(tan(theta / 2.));

  return eta;
}

double RawTowerGeomv1::get_phi() const
{
  return atan2(_center_y, _center_x);
}

void RawTowerGeomv1::identify(std::ostream& os) const
{
  os << "RawTowerGeomv1:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z() << std::endl;
}
