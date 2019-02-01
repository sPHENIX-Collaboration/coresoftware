#include "RawTowerGeomv3.h"

#include <cmath>
#include <iostream>

using namespace std;

RawTowerGeomv3::RawTowerGeomv3()
  : _towerid(~0)
  , _center_x(0)
  , _center_y(0)
  , _center_z(0)
  , _size_x(0)
  , _size_y(0)
  , _size_z(0)
  , _tower_type(-1)
{
}

RawTowerGeomv3::RawTowerGeomv3(RawTowerDefs::keytype id)
  : _towerid(id)
  , _center_x(0)
  , _center_y(0)
  , _center_z(0)
  , _size_x(0)
  , _size_y(0)
  , _size_z(0)
  , _tower_type(-1)
{
}

double RawTowerGeomv3::get_center_radius() const
{
  return sqrt(_center_x * _center_x +
              _center_y * _center_y);
}

double RawTowerGeomv3::get_eta() const
{
  double eta;
  double radius;
  double theta;
  radius = sqrt(_center_x * _center_x + _center_y * _center_y);
  theta = atan2(radius, _center_z);
  eta = -log(tan(theta / 2.));
  return eta;
}

double RawTowerGeomv3::get_phi() const
{
  return atan2(_center_y, _center_x);
}

void RawTowerGeomv3::identify(std::ostream& os) const
{
  os << "RawTowerGeomv3:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z()
     << "\n           dx: " << get_size_x() << " dy: " << get_size_y() << " dz: " << get_size_z()
     << "\n           tower_type = " << _tower_type << std::endl;
}
