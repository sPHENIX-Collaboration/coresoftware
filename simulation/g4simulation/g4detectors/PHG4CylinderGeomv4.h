// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOMV4_H
#define G4DETECTORS_PHG4CYLINDERGEOMV4_H

#include "PHG4CylinderGeom.h"

#include <iostream>  // for cout, ostream

class PHG4CylinderGeomv4 : public PHG4CylinderGeom
{
 public:
  PHG4CylinderGeomv4();
  PHG4CylinderGeomv4(const int lnsensors,
                     const int lnz,
                     const int nspc,
                     int nsc,
                     const int nstag,
                     const double lr,
                     const double rs,
                     const double szs,
                     const double sps,
                     const double sxo,
                     double syo,
                     const double szsp,
                     const double sys,
                     const double tck,
                     const double st)
    : N_sensors_in_layer(lnsensors)
    , layer(-1)
    , layer_radius(lr)
    , radius_stagger(rs)
    , layer_NZ(lnz)
    , segment_z_step(szs)
    , segment_phi_step(sps)
    , sensor_x_offset(sxo)
    , sensor_y_offset(syo)
    , N_strip_columns(nsc)
    , N_strips_per_column(nspc)
    , N_staggers(nstag)
    , strip_z_spacing(szsp)
    , strip_y_spacing(sys)
    , thickness(tck)
    , strip_tilt(st)
  {
  }

  ~PHG4CylinderGeomv4() override {}

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  void set_layer(const int i) override { layer = i; }
  int get_layer() const override { return layer; }
  double get_radius() const override { return layer_radius; }

  void find_segment_center(const int segment_z_bin, const int segment_phi_bin, double location[]) override;
  void find_strip_center(const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[]) override;

  double get_thickness() const override { return thickness; }
  double get_strip_y_spacing() const override { return strip_y_spacing; }
  double get_strip_z_spacing() const override { return strip_z_spacing; }
  double get_strip_tilt() const override { return strip_tilt; }
  int get_N_strip_columns() const override { return N_strip_columns; }
  int get_N_strips_per_column() const override { return N_strips_per_column; }
  int get_N_sensors_in_layer() const override { return N_sensors_in_layer; }

  // our own (not inherited from base class)
  double get_sensor_x_offset() const { return sensor_x_offset; }
  double get_sensor_y_offset() const { return sensor_y_offset; }

 protected:
  int N_sensors_in_layer;
  int layer;

  // finding the center of a sensor given ladder_segment_z and ladder_
  double layer_radius;
  double radius_stagger;
  int layer_NZ;
  double segment_z_step;
  double segment_phi_step;
  double sensor_x_offset;
  double sensor_y_offset;

  // navigation within a sensor
  //double strip_x_offset;
  // double strip_y_offset;
  int N_strip_columns;
  int N_strips_per_column;
  int N_staggers;
  double strip_z_spacing;
  double strip_y_spacing;
  double thickness;
  double strip_tilt;

  ClassDefOverride(PHG4CylinderGeomv4, 1)
};

#endif
