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

  _rot[0] = geom0.get_rotx();
  _rot[1] = geom0.get_roty();
  _rot[2] = geom0.get_rotz();

  for (int i = 0; i < _nVtx; i++)
  {
    _vertices[(i * _nDim) + 0] = geom0.get_vertex_x(i);
    _vertices[(i * _nDim) + 1] = geom0.get_vertex_y(i);
    _vertices[(i * _nDim) + 2] = geom0.get_vertex_z(i);
  }

  _center[0] = geom0.get_center_x();
  _center[1] = geom0.get_center_y();
  _center[2] = geom0.get_center_z();

  _center_int[0] = geom0.get_center_int_x();
  _center_int[1] = geom0.get_center_int_y();
  _center_int[2] = geom0.get_center_int_z();

  _center_ext[0] = geom0.get_center_ext_x();
  _center_ext[1] = geom0.get_center_ext_y();
  _center_ext[2] = geom0.get_center_ext_z();

  _center_low_eta[0] = geom0.get_center_low_eta_x();
  _center_low_eta[1] = geom0.get_center_low_eta_y();
  _center_low_eta[2] = geom0.get_center_low_eta_z();

  _center_high_eta[0] = geom0.get_center_high_eta_x();
  _center_high_eta[1] = geom0.get_center_high_eta_y();
  _center_high_eta[2] = geom0.get_center_high_eta_z();

  _center_low_phi[0] = geom0.get_center_low_phi_x();
  _center_low_phi[1] = geom0.get_center_low_phi_y();
  _center_low_phi[2] = geom0.get_center_low_phi_z();

  _center_high_phi[0] = geom0.get_center_high_phi_x();
  _center_high_phi[1] = geom0.get_center_high_phi_y();
  _center_high_phi[2] = geom0.get_center_high_phi_z();

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

  _center.fill(0.0);
  _center_int.fill(0.0);
  _center_ext.fill(0.0);
  _center_low_eta.fill(0.0);
  _center_high_eta.fill(0.0);
  _center_low_phi.fill(0.0);
  _center_high_phi.fill(0.0);
  
  for (int iDim = 0; iDim < _nDim; iDim++)
  {
    for (int iVtx = 0; iVtx < _nVtx; iVtx++)
    {
      _vertices[(iVtx * _nDim) + iDim] = vertices[(iVtx * _nDim) + iDim];
      
      _center[iDim] += _vertices[(iVtx * _nDim) + iDim];

      if (face_inside[iVtx])
      {
        _center_int[iDim] +=  _vertices[(iVtx * _nDim) + iDim];
      }
      if (face_outside[iVtx])
      {
        _center_ext[iDim] +=  _vertices[(iVtx * _nDim) + iDim];
      }
      if (face_low_eta[iVtx])
      {
        _center_low_eta[iDim] +=  _vertices[(iVtx * _nDim) + iDim];
      }
      if (face_high_eta[iVtx])
      {
        _center_high_eta[iDim] +=  _vertices[(iVtx * _nDim) + iDim];
      }
      if (face_low_phi[iVtx])
      {
        _center_low_phi[iDim] +=  _vertices[(iVtx * _nDim) + iDim];
      }
      if (face_high_phi[iVtx])
      {
        _center_high_phi[iDim] +=  _vertices[(iVtx * _nDim) + iDim];
      }
    }

    _center[iDim] /= (double) _nVtx;

    _center_int[iDim] /= (double) _nVtx / 2.;

    _center_ext[iDim] /= (double) _nVtx / 2.;
    
    _center_low_eta[iDim] /= (double) _nVtx / 2.;
    
    _center_high_eta[iDim] /= (double) _nVtx / 2.;
    
    _center_low_phi[iDim] /= (double) _nVtx / 2.;
    
    _center_high_phi[iDim] /= (double) _nVtx / 2.;
  }

  return;
}

double RawTowerGeomv5::get_vertex_x(int i) const
{
  if (i < 0 || i >= _nVtx)
  {
    std::cerr << "RawTowerGeomv5::get_vertex_x(): Vertex index " << i << " out of bounds. Should be between 0 and " << _nVtx << std::endl;
    exit(1);
  }
  return _vertices[i * _nDim + 0];
}

double RawTowerGeomv5::get_vertex_y(int i) const
{
  if (i < 0 || i >= _nVtx)
  {
    std::cerr << "RawTowerGeomv5::get_vertex_y(): Vertex index " << i << " out of bounds. Should be between 0 and " << _nVtx << std::endl;
    exit(1);
  }
  return _vertices[i * _nDim + 1];
}

double RawTowerGeomv5::get_vertex_z(int i) const
{
  if (i < 0 || i >= _nVtx)
  {
    std::cerr << "RawTowerGeomv5::get_vertex_z(): Vertex index " << i << " out of bounds. Should be between 0 and " << _nVtx << std::endl;
    exit(1);
  }
  return _vertices[i * _nDim + 2];
}

double RawTowerGeomv5::get_center_radius() const
{
  return std::sqrt(_center[0] * _center[0] +
                   _center[1] * _center[1]);
}

double RawTowerGeomv5::get_theta() const
{
  double radius = std::sqrt(_center[0] * _center[0] + _center[1] * _center[1]);
  double theta = std::atan2(radius, _center[2]);
  return theta;
}

double RawTowerGeomv5::get_eta() const
{
  double theta = get_theta();
  double eta = -std::log(std::tan(theta / 2.));

  return eta;
}

double RawTowerGeomv5::get_phi() const
{
  return std::atan2(_center[1], _center[0]);
}

void RawTowerGeomv5::identify(std::ostream& os) const
{
  os << "RawTowerGeomv5:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z() << std::endl;
}
