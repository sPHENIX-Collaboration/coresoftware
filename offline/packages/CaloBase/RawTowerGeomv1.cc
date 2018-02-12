#include "RawTowerGeomv1.h"

#include <cmath>
#include <iostream>
#include <limits>

using namespace std;

RawTowerGeomv1::RawTowerGeomv1()
  : RawTowerGeomv1(0)
{
}

RawTowerGeomv1::RawTowerGeomv1(RawTowerDefs::keytype id)
  : _towerid(id)
  , _center_x(numeric_limits<double>::signaling_NaN())
  , _center_y(numeric_limits<double>::signaling_NaN())
  , _center_z(numeric_limits<double>::signaling_NaN())
  , _tower_lenth_vec_x(numeric_limits<double>::signaling_NaN())
  , _tower_lenth_vec_y(numeric_limits<double>::signaling_NaN())
  , _tower_lenth_vec_z(numeric_limits<double>::signaling_NaN())
{
}

double RawTowerGeomv1::get_center_radius() const
{
  return sqrt(_center_x * _center_x +
              _center_y * _center_y);
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
     << "\n           length x: " << get_tower_lenth_vec_x() << " length y: " << get_tower_lenth_vec_y() << " length z: " << get_tower_lenth_vec_z() << std::endl;
}

double RawTowerGeomv1::get_front_x() const { return get_tower_lenth_vec_x() - 0.5 * get_tower_lenth_vec_x(); }
double RawTowerGeomv1::get_front_y() const { return get_tower_lenth_vec_y() - 0.5 * get_tower_lenth_vec_y(); }
double RawTowerGeomv1::get_front_z() const { return get_tower_lenth_vec_z() - 0.5 * get_tower_lenth_vec_z(); }
