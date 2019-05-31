#ifndef G4DETECTORS_PHG4BLOCKGEOMV1_H
#define G4DETECTORS_PHG4BLOCKGEOMV1_H

#include "PHG4BlockGeom.h"

#include <iostream>         // for cout, ostream

class PHG4BlockGeomv1: public PHG4BlockGeom
{
 public:
  PHG4BlockGeomv1();
  PHG4BlockGeomv1( const int layer,
                   const double sizex, const double sizey, const double sizez, 
                   const double centerx, const double centery, const double centerz,
                   const double zrot );

  virtual ~PHG4BlockGeomv1() {}

  void identify(std::ostream& os = std::cout) const;

  int get_layer() const {return _layer;}
  double get_width() const {return _size[0];}
  double get_thickness() const {return _size[1];}
  double get_length() const {return _size[2];}
  double get_center_x() const {return _center[0];}
  double get_center_y() const {return _center[1];}
  double get_center_z() const {return _center[2];}
  double get_z_rot() const {return _rotation_z;}

  double get_size_x() const {return _size[0];}
  double get_size_y() const {return _size[1];}
  double get_size_z() const {return _size[2];}

  double get_rot_matrix(const int i, const int j) const {return _rot_matrix[i][j];}

  void set_layer(const int i) {_layer = i;}

  // size in local coordinates
  void set_size(const double sizex, const double sizey, const double sizez);

  void set_center(const double centerx, const double centery, const double centerz);
  void set_z_rot(const double zrot) {_build_rot_matrix(); _rotation_z = zrot;}

  void convert_local_to_global(double, double, double,
                               double &, double &, double &) const;
  void convert_global_x_to_local(double, double, double, 
                                 double &, double &, double &) const;

 protected:
  int _layer;
  double _size[3];
  double _center[3];
  double _rotation_z;

  void _build_rot_matrix();
  double _rot_matrix[3][3];  // global -> local coordinates rotation matrix

  ClassDef(PHG4BlockGeomv1,1)
};

#endif
