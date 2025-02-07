#ifndef CALOBASE_RAWTOWERGEO_H
#define CALOBASE_RAWTOWERGEO_H

#include "RawTowerDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <iostream>
#include <limits>

class RawTowerGeom : public PHObject
{
 public:
  ~RawTowerGeom() override = default;

  //RawTowerGeom(const RawTowerGeom& /*geom*/){}

  void identify(std::ostream& os = std::cout) const override;

  virtual void set_id(RawTowerDefs::keytype /*key*/) { PHOOL_VIRTUAL_WARN("set_id()"); }

  virtual RawTowerDefs::keytype get_id() const
  {
    PHOOL_VIRTUAL_WARN("get_id()");
    return 0;
  }

  virtual int get_bineta() const
  {
    PHOOL_VIRTUAL_WARN("get_ieta()");
    return -1;
  }

  virtual int get_binphi() const
  {
    PHOOL_VIRTUAL_WARN("get_iphi()");
    return -1;
  }

  virtual int get_binl() const
  {
    PHOOL_VIRTUAL_WARN("get_binl()");
    return -1;
  }

  virtual int get_column() const
  {
    PHOOL_VIRTUAL_WARN("get_column()");
    return -1;
  }

  virtual int get_row() const
  {
    PHOOL_VIRTUAL_WARN("get_row()");
    return -1;
  }

  virtual void set_center_x(double /*x*/)
  {
    PHOOL_VIRTUAL_WARN("set_center_x()");
    return;
  }

  virtual void set_center_y(double /*y*/)
  {
    PHOOL_VIRTUAL_WARN("set_center_y()");
    return;
  }

  virtual void set_center_z(double /*z*/)
  {
    PHOOL_VIRTUAL_WARN("set_center_z()");
    return;
  }

  virtual void set_vertices(const std::vector<double>& /*vertices*/)
  {
    PHOOL_VIRTUAL_WARN("set_vertices()");
    return;
  }

  virtual void set_size_x(double /*dx*/)
  {
    PHOOL_VIRTUAL_WARN("set_size_x()");
    return;
  }

  virtual void set_size_y(double /*dy*/)
  {
    PHOOL_VIRTUAL_WARN("set_size_y()");
    return;
  }

  virtual void set_size_z(double /*dz*/)
  {
    PHOOL_VIRTUAL_WARN("set_size_z()");
    return;
  }

  virtual double get_center_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_int_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_int_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_int_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_int_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_int_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_int_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_ext_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_ext_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_ext_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_ext_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_ext_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_ext_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_low_eta_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_low_eta_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_low_eta_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_low_eta_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_low_eta_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_low_eta_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_high_eta_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_high_eta_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_high_eta_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_high_eta_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_high_eta_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_high_eta_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_low_phi_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_low_phi_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_low_phi_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_low_phi_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_low_phi_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_low_phi_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_high_phi_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_high_phi_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_high_phi_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_high_phi_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_high_phi_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_high_phi_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_vertex_x(int /*i*/) const
  {
    PHOOL_VIRTUAL_WARN("get_vertex_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_vertex_y(int /*i*/) const
  {
    PHOOL_VIRTUAL_WARN("get_vertex_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_vertex_z(int /*i*/) const
  {
    PHOOL_VIRTUAL_WARN("get_vertex_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_size_x() const
  {
    PHOOL_VIRTUAL_WARN("get_size_x()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_size_y() const
  {
    PHOOL_VIRTUAL_WARN("get_size_y()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_size_z() const
  {
    PHOOL_VIRTUAL_WARN("get_size_z()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_volume() const
  {
    PHOOL_VIRTUAL_WARN("get_volume()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_center_radius() const
  {
    PHOOL_VIRTUAL_WARN("get_center_radius()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_eta() const
  {
    PHOOL_VIRTUAL_WARN("get_eta()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_theta() const
  {
    PHOOL_VIRTUAL_WARN("get_theta()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual double get_phi() const
  {
    PHOOL_VIRTUAL_WARN("get_phi()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual void set_tower_type(int /*tt*/)
  {
    PHOOL_VIRTUAL_WARN("set_tower_type()");
    return;
  }
  virtual int get_tower_type() const
  {
    PHOOL_VIRTUAL_WARN("get_tower_type()");
    return -1;
  }

  virtual double get_rotx() const
  {
    PHOOL_VIRTUAL_WARN("get_rotx()");
    return std::numeric_limits<float>::quiet_NaN();
  }
  virtual double get_roty() const
  {
    PHOOL_VIRTUAL_WARN("get_roty()");
    return std::numeric_limits<float>::quiet_NaN();
  }
  virtual double get_rotz() const
  {
    PHOOL_VIRTUAL_WARN("get_rotz()");
    return std::numeric_limits<float>::quiet_NaN();
  }

  virtual void set_rotx(double /*rotx*/)
  {
    PHOOL_VIRTUAL_WARN("set_rotx()");
    return;
  }
  virtual void set_roty(double /*roty*/)
  {
    PHOOL_VIRTUAL_WARN("set_roty()");
    return;
  }
  virtual void set_rotz(double /*rotz*/)
  {
    PHOOL_VIRTUAL_WARN("set_rotz()");
    return;
  }

 protected:
  RawTowerGeom() = default;

  ClassDefOverride(RawTowerGeom, 2)
};

#endif
