#include "RawTowerGeomv1.h"

#include <iostream>
#include <cmath>

using namespace std;

RawTowerGeomv1::RawTowerGeomv1() :
  _towerid(~0),
  _center_x(0),
  _center_y(0),
  _center_z(0)
{}

RawTowerGeomv1::RawTowerGeomv1(RawTowerDefs::keytype id) :
  _towerid(id),
  _center_x(0),
  _center_y(0),
  _center_z(0)
{}

double RawTowerGeomv1::get_center_radius() const
{
  return sqrt( _center_x * _center_x +
	       _center_y * _center_y );
}

double RawTowerGeomv1::get_eta() const
{


  double eta;
  double radius;
  double theta;
  radius = sqrt(_center_x * _center_x + _center_y * _center_y);
  theta = atan2(radius, _center_z);
  eta = -log(tan(theta / 2.));

  return eta;
}

double RawTowerGeomv1::get_phi() const
{
  return atan2(_center_y, _center_x);
}


void RawTowerGeomv1::identify(std::ostream& os) const
{
  os << "RawTowerGeomv1:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z()
	    << "\n           dx: " << get_size_x() << " dy: " << get_size_y() << " dz: " << get_size_z() << std::endl;
}
