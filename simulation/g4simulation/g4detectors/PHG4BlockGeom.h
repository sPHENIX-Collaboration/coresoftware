// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BLOCKGEOM_H
#define G4DETECTORS_PHG4BLOCKGEOM_H

#include <phool/PHObject.h>

#include <phool/phool.h>

#include <limits>
#include <iostream>  // for cout, ostream

class PHG4BlockGeom : public PHObject
{
 public:
  ~PHG4BlockGeom() = default;

  // from PHObject
  void identify(std::ostream &os = std::cout) const override;

  virtual int get_layer() const
  {
    PHOOL_VIRTUAL_WARN("get_layer()");
    return -99999;
  }
  virtual double get_size_x() const
  {
    PHOOL_VIRTUAL_WARN("get_size_x()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_size_y() const
  {
    PHOOL_VIRTUAL_WARN("get_size_y()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_size_z() const
  {
    PHOOL_VIRTUAL_WARN("get_size_z()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_center_x() const
  {
    PHOOL_VIRTUAL_WARN("get_place_x()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_center_y() const
  {
    PHOOL_VIRTUAL_WARN("get_place_y()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_center_z() const
  {
    PHOOL_VIRTUAL_WARN("get_place_z()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_z_rot() const
  {
    PHOOL_VIRTUAL_WARN("get_z_rot()");
    return std::numeric_limits<double>::quiet_NaN();
  }

  virtual double get_width() const
  {
    PHOOL_VIRTUAL_WARN("get_width()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_length() const
  {
    PHOOL_VIRTUAL_WARN("get_length()");
    return std::numeric_limits<double>::quiet_NaN();
  }
  virtual double get_thickness() const
  {
    PHOOL_VIRTUAL_WARN("get_thickness()");
    return std::numeric_limits<double>::quiet_NaN();
  }

  virtual double get_rot_matrix(const int, const int) const
  {
    PHOOL_VIRTUAL_WARN("get_rot_matrix(const int, const int)");
    return std::numeric_limits<double>::quiet_NaN();
  }

  virtual void set_layer(const int) { PHOOL_VIRTUAL_WARN("set_layer(const int)"); }
  virtual void set_size(const double /*sizex*/, const double /*sizey*/, const double /*sizez*/)
  {
    PHOOL_VIRTUAL_WARN("set_size(const double, const double, const double)");
  }
  virtual void set_place(const double /*placex*/, const double /*placey*/, const double /*placez*/)
  {
    PHOOL_VIRTUAL_WARN("set_place(const double, const double, const double)");
  }
  virtual void set_z_rot(const double) { PHOOL_VIRTUAL_WARN("set_z_rot(const double)"); }

  virtual void convert_local_to_global(const double, const double, const double,
                                       double &, double &, double &) const
  {
    PHOOL_VIRTUAL_WARN("convert_local_to_global(const double, const double, const double, double &, double &, double &)");
  }
  virtual void convert_global_to_local(const double, const double, const double,
                                       double &, double &, double &) const
  {
    PHOOL_VIRTUAL_WARN("convert_global_to_local(const double, const double, const double, double &, double &, double &)");
  }

 protected:
  PHG4BlockGeom() = default;

  ClassDefOverride(PHG4BlockGeom, 1)
};

#endif
