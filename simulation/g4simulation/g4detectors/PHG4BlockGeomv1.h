// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKGEOMV1_H
#define G4DETECTORS_PHG4BLOCKGEOMV1_H

#include "PHG4BlockGeom.h"

#include <iostream>  // for cout, ostream

class PHG4BlockGeomv1 : public PHG4BlockGeom
{
 public:
  PHG4BlockGeomv1();
  PHG4BlockGeomv1(const int layer,
                  const double sizex, const double sizey, const double sizez,
                  const double centerx, const double centery, const double centerz,
                  const double zrot);

  ~PHG4BlockGeomv1() override {}

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  int get_layer() const override { return _layer; }
  double get_width() const override { return _size[0]; }
  double get_thickness() const override { return _size[1]; }
  double get_length() const override { return _size[2]; }
  double get_center_x() const override { return _center[0]; }
  double get_center_y() const override { return _center[1]; }
  double get_center_z() const override { return _center[2]; }
  double get_z_rot() const override { return _rotation_z; }

  double get_size_x() const override { return _size[0]; }
  double get_size_y() const override { return _size[1]; }
  double get_size_z() const override { return _size[2]; }

  double get_rot_matrix(const int i, const int j) const override { return _rot_matrix[i][j]; }

  void set_layer(const int i) override { _layer = i; }

  // size in local coordinates
  void set_size(const double sizex, const double sizey, const double sizez) override;

  void set_z_rot(const double zrot) override
  {
    _build_rot_matrix();
    _rotation_z = zrot;
  }

  void convert_local_to_global(double, double, double,
                               double &, double &, double &) const override;

  // our own (not inherited)
  void set_center(const double centerx, const double centery, const double centerz);
  void convert_global_x_to_local(double, double, double,
                                 double &, double &, double &) const;

 protected:
  int _layer;
  double _size[3];
  double _center[3];
  double _rotation_z;

  void _build_rot_matrix();
  double _rot_matrix[3][3];  // global -> local coordinates rotation matrix

  ClassDefOverride(PHG4BlockGeomv1, 1)
};

#endif
