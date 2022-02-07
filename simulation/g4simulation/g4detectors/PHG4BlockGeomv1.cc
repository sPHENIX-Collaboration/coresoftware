#include "PHG4BlockGeomv1.h"

#include <algorithm>
#include <cmath>

using namespace std;

PHG4BlockGeomv1::PHG4BlockGeomv1()
  : PHG4BlockGeom()
  , _layer(-1)
  , _rotation_z(NAN)
{
  const double filldval = 1.;
  fill(_size, _size + sizeof(_size) / sizeof(double), NAN);
  fill(_center, _center + sizeof(_center) / sizeof(double), NAN);
  fill(&_rot_matrix[0][0], &_rot_matrix[0][0] + sizeof(_rot_matrix) / sizeof(double), filldval);
}

PHG4BlockGeomv1::PHG4BlockGeomv1(const int layer,
                                 const double sizex, const double sizey, const double sizez,
                                 const double centerx, const double centery, const double centerz,
                                 const double zrot)
  : PHG4BlockGeom()
  , _layer(layer)
  , _rotation_z(zrot)
{
  _size[0] = sizex;
  _size[1] = sizey;
  _size[2] = sizez;
  _center[0] = centerx;
  _center[1] = centery;
  _center[2] = centerz;

  _build_rot_matrix();
}

void PHG4BlockGeomv1::identify(std::ostream &os) const
{
  os << "PHG4BlockGeomv1: layer: " << _layer
     << ", rotation in z: " << _rotation_z
     << ", size: (" << _size[0] << ", " << _size[1] << ", " << _size[2] << ")"
     << ", center: (" << _center[0] << ", " << _center[1] << ", " << _center[2] << ")"
     << endl;
  return;
}

void PHG4BlockGeomv1::set_size(const double sizex, const double sizey, const double sizez)
{
  _size[0] = sizex;
  _size[1] = sizey;
  _size[2] = sizez;
  return;
}

void PHG4BlockGeomv1::set_center(const double centerx, const double centery, const double centerz)
{
  _center[0] = centerx;
  _center[1] = centery;
  _center[2] = centerz;
  return;
}

void PHG4BlockGeomv1::convert_local_to_global(double lx, double ly, double lz,
                                              double &gx, double &gy, double &gz) const
{
  // gvec = R^T lvec + offset
  gx = _rot_matrix[0][0] * lx + _rot_matrix[1][0] * ly + _rot_matrix[2][0] * lz;
  gy = _rot_matrix[0][1] * lx + _rot_matrix[1][1] * ly + _rot_matrix[2][1] * lz;
  gz = _rot_matrix[0][2] * lx + _rot_matrix[1][2] * ly + _rot_matrix[2][2] * lz;
  gx += _center[0];
  gy += _center[1];
  gz += _center[2];
  return;
}

void PHG4BlockGeomv1::convert_global_x_to_local(double gx, double gy, double gz,
                                                double &lx, double &ly, double &lz) const
{
  // lvec = R (gvec - offset)
  gx -= _center[0];
  gy -= _center[1];
  gz -= _center[2];
  lx = _rot_matrix[0][0] * gx + _rot_matrix[0][1] * gy + _rot_matrix[0][2] * gz;
  ly = _rot_matrix[1][0] * gx + _rot_matrix[1][1] * gy + _rot_matrix[1][2] * gz;
  lz = _rot_matrix[2][0] * gx + _rot_matrix[2][1] * gy + _rot_matrix[2][2] * gz;
  return;
}

void PHG4BlockGeomv1::_build_rot_matrix()
{
  _rot_matrix[0][0] = cos(_rotation_z);
  _rot_matrix[0][1] = sin(_rotation_z);
  _rot_matrix[1][0] = -sin(_rotation_z);
  _rot_matrix[1][1] = cos(_rotation_z);
  _rot_matrix[2][2] = 1.;
}
