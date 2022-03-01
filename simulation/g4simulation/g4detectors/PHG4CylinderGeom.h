// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CYLINDERGEOM_H
#define G4DETECTORS_PHG4CYLINDERGEOM_H

#include <phool/PHObject.h>

#include <phool/phool.h>

#include <cmath>
#include <iostream>  // for cout, ostream

class PHParameters;

class PHG4CylinderGeom : public PHObject
{
 public:
  ~PHG4CylinderGeom() override {}

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  virtual int get_layer() const
  {
    PHOOL_VIRTUAL_WARN("get_layer()");
    return -99999;
  }
  virtual double get_radius() const
  {
    PHOOL_VIRTUAL_WARN("get_radius()");
    return NAN;
  }
  virtual double get_thickness() const
  {
    PHOOL_VIRTUAL_WARN("get_thickness()");
    return NAN;
  }
  virtual double get_zmin() const
  {
    PHOOL_VIRTUAL_WARN("get_zmin()");
    return NAN;
  }
  virtual double get_zmax() const
  {
    PHOOL_VIRTUAL_WARN("get_zmax()");
    return NAN;
  }
  virtual int get_nscint() const
  {
    PHOOL_VIRTUAL_WARN("get_nscint()");
    return -99999;
  }
  virtual double get_tiltangle() const
  {
    PHOOL_VIRTUAL_WARN("get_tiltangle()");
    return NAN;
  }
  virtual double get_phi_slat_zero() const
  {
    PHOOL_VIRTUAL_WARN("get_phi_slat_zero()");
    return NAN;
  }

  virtual void set_layer(const int) { PHOOL_VIRTUAL_WARN("set_layer(const int)"); }
  virtual void set_radius(const double) { PHOOL_VIRTUAL_WARN("set_radius(const double)"); }
  virtual void set_thickness(const double) { PHOOL_VIRTUAL_WARN("set_thickness(const double)"); }
  virtual void set_zmin(const double) { PHOOL_VIRTUAL_WARN("set_zmin(const double)"); }
  virtual void set_zmax(const double) { PHOOL_VIRTUAL_WARN("set_zmax(const double)"); }
  virtual void set_nscint(const int)
  {
    PHOOL_VIRTUAL_WARN("set_nscint(const int)");
    return;
  }
  virtual void set_tiltangle(const double)
  {
    PHOOL_VIRTUAL_WARN("set_tiltangle(const double)");
    return;
  }
  virtual void set_phi_slat_zero(const double)
  {
    PHOOL_VIRTUAL_WARN("set_phi_slat_zero (const double)");
    return;
  }

  virtual void find_segment_center(const int /*segment_z_bin*/, const int /*segment_phi_bin*/, double /*location*/[])
  {
    PHOOL_VIRTUAL_WARN("find_sensor_center");
    return;
  }
  virtual void find_strip_center(const int /*segment_z_bin*/, const int /*segment_phi_bin*/, const int /*strip_column*/, const int /*strip_index*/, double /*location*/[])
  {
    PHOOL_VIRTUAL_WARN("find_strip_center");
    return;
  }
  virtual void find_strip_index_values(const int /*segment_z_bin*/, const double /*ypos*/, const double /*zpos*/, int & /*strip_y_index*/, int & /*strip_z_index*/)
  {
    PHOOL_VIRTUAL_WARN("find_strip_index_values");
    return;
  }
  virtual void find_strip_center_local_coords(const int /*segment_z_bin*/, const int /*strip_y_index*/, const int /*strip_z_index*/, double /*location*/[])
  {
    PHOOL_VIRTUAL_WARN("find_strip_center_localcoords");
    return;
  }

  virtual double get_strip_y_spacing() const
  {
    PHOOL_VIRTUAL_WARN("get_strip_y_spacing");
    return NAN;
  }
  virtual double get_strip_z_spacing() const
  {
    PHOOL_VIRTUAL_WARN("get_strip_z_spacing");
    return NAN;
  }
  virtual double get_strip_tilt() const
  {
    PHOOL_VIRTUAL_WARN("get_strip_tilt");
    return NAN;
  }

  virtual int get_N_strip_columns() const
  {
    PHOOL_VIRTUAL_WARN("get_N_strip_columns");
    return -9999;
  }
  virtual int get_N_strips_per_column() const
  {
    PHOOL_VIRTUAL_WARN("get_N_strips_per_column");
    return -9999;
  }
  virtual int get_N_sensors_in_layer() const
  {
    PHOOL_VIRTUAL_WARN("get_N_sensors_in_layer");
    return -9999;
  }

  virtual double get_pixel_z() const
  {
    PHOOL_VIRTUAL_WARN("get_pixel_z");
    return NAN;
  }
  virtual double get_pixel_x() const
  {
    PHOOL_VIRTUAL_WARN("get_pixel_x");
    return NAN;
  }
  virtual double get_pixel_thickness() const
  {
    PHOOL_VIRTUAL_WARN("get_pixel_thickness");
    return NAN;
  }

  //! load parameters from PHParameters, which interface to Database/XML/ROOT files
  virtual void ImportParameters(const PHParameters & /*param*/) { return; }

 protected:
  PHG4CylinderGeom() {}

  ClassDefOverride(PHG4CylinderGeom, 1)
};

#endif
