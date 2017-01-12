#include "RawTowerGeomv2.h"

#include <iostream>
#include <cmath>

using namespace std;

RawTowerGeomv2::RawTowerGeomv2() :
  _towerid(~0),
  _center_x(0),
  _center_y(0),
  _center_z(0),
  _size_x(0),
  _size_y(0),
  _size_z(0)
{}
RawTowerGeomv2::RawTowerGeomv2(RawTowerDefs::keytype id) :
  _towerid(id),
  _center_x(0),
  _center_y(0),
  _center_z(0),
  _size_x(0),
  _size_y(0),
  _size_z(0)
{}

double RawTowerGeomv2::get_center_radius() const
{
  return sqrt( _center_x * _center_x +
	       _center_y * _center_y );
}

double RawTowerGeomv2::get_eta() const
{
  double eta;
  double radius;
  double theta;
  radius = sqrt(_center_x * _center_x + _center_y * _center_y);
  theta = atan2(radius, _center_z);
  eta = -log(tan(theta / 2.));

  return eta;
}

double RawTowerGeomv2::get_phi() const
{
  return atan2(_center_y, _center_x);
}

void RawTowerGeomv2::identify(std::ostream& os) const
{
  os << "RawTowerGeomv2:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z()
	    << "\n           dx: " << get_size_x() << " dy: " << get_size_y() << " dz: " << get_size_z() << std::endl;
}
