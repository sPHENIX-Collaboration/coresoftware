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
  
  for (int i = 0; i < _nVtx; i++)
  {
    _vertices_x[i] = geom0.get_vertex_x(i);
    _vertices_y[i] = geom0.get_vertex_y(i);
    _vertices_z[i] = geom0.get_vertex_z(i);
  }

  _center_x = geom0.get_center_x();
  _center_y = geom0.get_center_y();
  _center_z = geom0.get_center_z();

  _center_int_x = geom0.get_center_int_x();
  _center_int_y = geom0.get_center_int_y();
  _center_int_z = geom0.get_center_int_z();
  
  _center_ext_x = geom0.get_center_ext_x();
  _center_ext_y = geom0.get_center_ext_y();
  _center_ext_z = geom0.get_center_ext_z();

  _center_low_eta_x = geom0.get_center_low_eta_x();
  _center_low_eta_y = geom0.get_center_low_eta_y();
  _center_low_eta_z = geom0.get_center_low_eta_z();

  _center_high_eta_x = geom0.get_center_high_eta_x();
  _center_high_eta_y = geom0.get_center_high_eta_y();
  _center_high_eta_z = geom0.get_center_high_eta_z();

  _center_low_phi_x = geom0.get_center_low_phi_x();
  _center_low_phi_y = geom0.get_center_low_phi_y();
  _center_low_phi_z = geom0.get_center_low_phi_z();

  _center_high_phi_x = geom0.get_center_high_phi_x();
  _center_high_phi_y = geom0.get_center_high_phi_y();
  _center_high_phi_z = geom0.get_center_high_phi_z();
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

  bool face_outside[_nVtx] = {false, false, false, false, true, true, true, true};
  bool face_inside[_nVtx] = {true, true, true, true, false, false, false, false};
  bool face_low_eta[_nVtx] = {false, false, true, true, false, false, true, true};
  bool face_high_eta[_nVtx] = {true, true, false, false, true, true, false, false};
  bool face_low_phi[_nVtx] = {false, true, true, false, false, true, true, false};
  bool face_high_phi[_nVtx] = {true, false, false, true, true, false, false, true};

  _center_x = 0;
  _center_y = 0;
  _center_z = 0;
  _center_int_x = 0;
  _center_int_y = 0;
  _center_int_z = 0;
  _center_ext_x = 0;
  _center_ext_y = 0;
  _center_ext_z = 0;
  _center_low_eta_x = 0;
  _center_low_eta_y = 0;
  _center_low_eta_z = 0;
  _center_high_eta_x = 0;
  _center_high_eta_y = 0;
  _center_high_eta_z = 0;
  _center_low_phi_x = 0;
  _center_low_phi_y = 0;
  _center_low_phi_z = 0;
  _center_high_phi_x = 0;
  _center_high_phi_y = 0;
  _center_high_phi_z = 0;
  for (int i = 0; i < _nVtx; i++)
  {
    _vertices_x[i] = vertices[i*3+0];
    _vertices_y[i] = vertices[i*3+1];
    _vertices_z[i] = vertices[i*3+2];

    _center_x += _vertices_x[i];
    _center_y += _vertices_y[i];
    _center_z += _vertices_z[i];

    if (face_inside[i])
    {
      _center_int_x +=  _vertices_x[i];
      _center_int_y +=  _vertices_y[i];
      _center_int_z +=  _vertices_z[i];
    }
    if (face_outside[i])
    {
      _center_ext_x +=  _vertices_x[i];
      _center_ext_y +=  _vertices_y[i];
      _center_ext_z +=  _vertices_z[i];
    }
    if (face_low_eta[i])
    {
      _center_low_eta_x +=  _vertices_x[i];
      _center_low_eta_y +=  _vertices_y[i];
      _center_low_eta_z +=  _vertices_z[i];
    }
    if (face_high_eta[i])
    {
      _center_high_eta_x +=  _vertices_x[i];
      _center_high_eta_y +=  _vertices_y[i];
      _center_high_eta_z +=  _vertices_z[i];
    }
    if (face_low_phi[i])
    {
      _center_low_phi_x +=  _vertices_x[i];
      _center_low_phi_y +=  _vertices_y[i];
      _center_low_phi_z +=  _vertices_z[i];
    }
    if (face_high_phi[i])
    {
      _center_high_phi_x +=  _vertices_x[i];
      _center_high_phi_y +=  _vertices_y[i];
      _center_high_phi_z +=  _vertices_z[i];
    }
  }

  _center_x /= (double) _nVtx;
  _center_y /= (double) _nVtx;
  _center_z /= (double) _nVtx;

  _center_int_x /= (double) _nVtx / 2.;
  _center_int_y /= (double) _nVtx / 2.;
  _center_int_z /= (double) _nVtx / 2.;

  _center_ext_x /= (double) _nVtx / 2.;
  _center_ext_y /= (double) _nVtx / 2.;
  _center_ext_z /= (double) _nVtx / 2.;
  
  _center_low_eta_x /= (double) _nVtx / 2.;
  _center_low_eta_y /= (double) _nVtx / 2.;
  _center_low_eta_z /= (double) _nVtx / 2.;

  _center_high_eta_x /= (double) _nVtx / 2.;
  _center_high_eta_y /= (double) _nVtx / 2.;
  _center_high_eta_z /= (double) _nVtx / 2.;

  _center_low_phi_x /= (double) _nVtx / 2.;
  _center_low_phi_y /= (double) _nVtx / 2.;
  _center_low_phi_z /= (double) _nVtx / 2.;

  _center_high_phi_x /= (double) _nVtx / 2.;
  _center_high_phi_y /= (double) _nVtx / 2.;
  _center_high_phi_z /= (double) _nVtx / 2.;

  return;
} 

double RawTowerGeomv5::get_center_radius() const
{
  return std::sqrt(_center_x * _center_x +
                   _center_y * _center_y);
}

double RawTowerGeomv5::get_theta() const
{
  double radius = sqrt(_center_x * _center_x + _center_y * _center_y);
  double theta = atan2(radius, _center_z);
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
  return atan2(_center_y, _center_x);
}

void RawTowerGeomv5::identify(std::ostream& os) const
{
  os << "RawTowerGeomv5:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z() << std::endl;
}
