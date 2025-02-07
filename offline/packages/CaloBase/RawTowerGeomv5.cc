#include "RawTowerGeomv5.h"

#include <cmath>
#include <iostream>

RawTowerGeomv5::RawTowerGeomv5(RawTowerDefs::keytype id)
  : _towerid(id)
{
}

RawTowerGeomv5::RawTowerGeomv5(const RawTowerGeom& geom0)
{
  _towerid = geom0.get_id();

  _rotx = geom0.get_rotx();
  _roty = geom0.get_roty();
  _rotz = geom0.get_rotz();

  _vertices.resize(_nVtx);
  for (int i = 0; i < _nVtx; i++)
  {
    _vertices[i].x = geom0.get_vertex_x(i);
    _vertices[i].y = geom0.get_vertex_y(i);
    _vertices[i].z = geom0.get_vertex_z(i);
  }

  _center.x = geom0.get_center_x();
  _center.y = geom0.get_center_y();
  _center.z = geom0.get_center_z();

  _center_int.x = geom0.get_center_int_x();
  _center_int.y = geom0.get_center_int_y();
  _center_int.z = geom0.get_center_int_z();

  _center_ext.x = geom0.get_center_ext_x();
  _center_ext.y = geom0.get_center_ext_y();
  _center_ext.z = geom0.get_center_ext_z();

  _center_low_eta.x = geom0.get_center_low_eta_x();
  _center_low_eta.y = geom0.get_center_low_eta_y();
  _center_low_eta.z = geom0.get_center_low_eta_z();

  _center_high_eta.x = geom0.get_center_high_eta_x();
  _center_high_eta.y = geom0.get_center_high_eta_y();
  _center_high_eta.z = geom0.get_center_high_eta_z();

  _center_low_phi.x = geom0.get_center_low_phi_x();
  _center_low_phi.y = geom0.get_center_low_phi_y();
  _center_low_phi.z = geom0.get_center_low_phi_z();

  _center_high_phi.x = geom0.get_center_high_phi_x();
  _center_high_phi.y = geom0.get_center_high_phi_y();
  _center_high_phi.z = geom0.get_center_high_phi_z();
}

void RawTowerGeomv5::set_vertices(const std::vector<double>& vertices)
{
  /*       y
           |                    7_______6
           |                   /|      /|
           |_______ z        8/_|_____/5|
          /                   | |     | |
         /                    |3|_____|_|2
        /                     | /     | /
     -x                      4|/______|/1
  */

  if ((int) vertices.size() != _nVtx * 3)
  {
    std::cerr
        << "RawTowerGeomv5::set_vertices - input " << vertices.size() << " vertices given. Expected 8 x 3." << std::endl;
    exit(1);
  }

  const bool face_outside[_nVtx] = {false, false, false, false, true, true, true, true};
  const bool face_inside[_nVtx] = {true, true, true, true, false, false, false, false};
  const bool face_low_eta[_nVtx] = {false, false, true, true, false, false, true, true};
  const bool face_high_eta[_nVtx] = {true, true, false, false, true, true, false, false};
  const bool face_low_phi[_nVtx] = {false, true, true, false, false, true, true, false};
  const bool face_high_phi[_nVtx] = {true, false, false, true, true, false, false, true};

  _center.x = 0;
  _center.y = 0;
  _center.z = 0;
  _center_int.x = 0;
  _center_int.y = 0;
  _center_int.z = 0;
  _center_ext.x = 0;
  _center_ext.y = 0;
  _center_ext.z = 0;
  _center_low_eta.x = 0;
  _center_low_eta.y = 0;
  _center_low_eta.z = 0;
  _center_high_eta.x = 0;
  _center_high_eta.y = 0;
  _center_high_eta.z = 0;
  _center_low_phi.x = 0;
  _center_low_phi.y = 0;
  _center_low_phi.z = 0;
  _center_high_phi.x = 0;
  _center_high_phi.y = 0;
  _center_high_phi.z = 0;
  _vertices.resize(_nVtx);
  for (int i = 0; i < _nVtx; i++)
  {
    _vertices[i].x = vertices[i * 3 + 0];
    _vertices[i].y = vertices[i * 3 + 1];
    _vertices[i].z = vertices[i * 3 + 2];

    _center.x += _vertices[i].x;
    _center.y += _vertices[i].y;
    _center.z += _vertices[i].z;

    if (face_inside[i])
    {
      _center_int.x += _vertices[i].x;
      _center_int.y += _vertices[i].y;
      _center_int.z += _vertices[i].z;
    }
    if (face_outside[i])
    {
      _center_ext.x += _vertices[i].x;
      _center_ext.y += _vertices[i].y;
      _center_ext.z += _vertices[i].z;
    }
    if (face_low_eta[i])
    {
      _center_low_eta.x += _vertices[i].x;
      _center_low_eta.y += _vertices[i].y;
      _center_low_eta.z += _vertices[i].z;
    }
    if (face_high_eta[i])
    {
      _center_high_eta.x += _vertices[i].x;
      _center_high_eta.y += _vertices[i].y;
      _center_high_eta.z += _vertices[i].z;
    }
    if (face_low_phi[i])
    {
      _center_low_phi.x += _vertices[i].x;
      _center_low_phi.y += _vertices[i].y;
      _center_low_phi.z += _vertices[i].z;
    }
    if (face_high_phi[i])
    {
      _center_high_phi.x += _vertices[i].x;
      _center_high_phi.y += _vertices[i].y;
      _center_high_phi.z += _vertices[i].z;
    }
  }

  _center.x /= (float) _nVtx;
  _center.y /= (float) _nVtx;
  _center.z /= (float) _nVtx;

  _center_int.x /= (float) _nVtx / 2.;
  _center_int.y /= (float) _nVtx / 2.;
  _center_int.z /= (float) _nVtx / 2.;

  _center_ext.x /= (float) _nVtx / 2.;
  _center_ext.y /= (float) _nVtx / 2.;
  _center_ext.z /= (float) _nVtx / 2.;

  _center_low_eta.x /= (float) _nVtx / 2.;
  _center_low_eta.y /= (float) _nVtx / 2.;
  _center_low_eta.z /= (float) _nVtx / 2.;

  _center_high_eta.x /= (float) _nVtx / 2.;
  _center_high_eta.y /= (float) _nVtx / 2.;
  _center_high_eta.z /= (float) _nVtx / 2.;

  _center_low_phi.x /= (float) _nVtx / 2.;
  _center_low_phi.y /= (float) _nVtx / 2.;
  _center_low_phi.z /= (float) _nVtx / 2.;

  _center_high_phi.x /= (float) _nVtx / 2.;
  _center_high_phi.y /= (float) _nVtx / 2.;
  _center_high_phi.z /= (float) _nVtx / 2.;

  return;
}

double RawTowerGeomv5::get_center_radius() const
{
  return std::sqrt(_center.x * _center.x +
                   _center.y * _center.y);
}

double RawTowerGeomv5::get_theta() const
{
  double radius = sqrt(_center.x * _center.x + _center.y * _center.y);
  double theta = atan2(radius, _center.z);
  return theta;
}

double RawTowerGeomv5::get_eta() const
{
  double theta = get_theta();
  double eta = -log(tan(theta / 2.));

  return eta;
}

double RawTowerGeomv5::get_phi() const
{
  return atan2(_center.y, _center.x);
}

void RawTowerGeomv5::identify(std::ostream& os) const
{
  os << "RawTowerGeomv5:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z() << std::endl;
}
