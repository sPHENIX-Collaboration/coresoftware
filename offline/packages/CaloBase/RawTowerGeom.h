#ifndef CALOBASE_RAWTOWERGEO_H
#define CALOBASE_RAWTOWERGEO_H

#include "RawTowerDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>

class RawTowerGeom : public PHObject
{
 public:
  ~RawTowerGeom() override {}

  void identify(std::ostream& os = std::cout) const override;

  virtual void set_id(RawTowerDefs::keytype) { PHOOL_VIRTUAL_WARN("set_id()"); }

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

  virtual void set_center_x(double)
  {
    PHOOL_VIRTUAL_WARN("set_center_x()");
    return;
  }

  virtual void set_center_y(double)
  {
    PHOOL_VIRTUAL_WARN("set_center_y()");
    return;
  }

  virtual void set_center_z(double)
  {
    PHOOL_VIRTUAL_WARN("set_center_z()");
    return;
  }

  virtual void set_size_x(double)
  {
    PHOOL_VIRTUAL_WARN("set_size_x()");
    return;
  }

  virtual void set_size_y(double)
  {
    PHOOL_VIRTUAL_WARN("set_size_y()");
    return;
  }

  virtual void set_size_z(double)
  {
    PHOOL_VIRTUAL_WARN("set_size_z()");
    return;
  }

  virtual double get_center_x() const
  {
    PHOOL_VIRTUAL_WARN("get_center_x()");
    return NAN;
  }

  virtual double get_center_y() const
  {
    PHOOL_VIRTUAL_WARN("get_center_y()");
    return NAN;
  }

  virtual double get_center_z() const
  {
    PHOOL_VIRTUAL_WARN("get_center_z()");
    return NAN;
  }

  virtual double get_size_x() const
  {
    PHOOL_VIRTUAL_WARN("get_size_x()");
    return NAN;
  }

  virtual double get_size_y() const
  {
    PHOOL_VIRTUAL_WARN("get_size_y()");
    return NAN;
  }

  virtual double get_size_z() const
  {
    PHOOL_VIRTUAL_WARN("get_size_z()");
    return NAN;
  }

  virtual double get_volume() const
  {
    PHOOL_VIRTUAL_WARN("get_volume()");
    return NAN;
  }

  virtual double get_center_radius() const
  {
    PHOOL_VIRTUAL_WARN("get_center_radius()");
    return NAN;
  }

  virtual double get_eta() const
  {
    PHOOL_VIRTUAL_WARN("get_eta()");
    return NAN;
  }

  virtual double get_theta() const
  {
    PHOOL_VIRTUAL_WARN("get_theta()");
    return NAN;
  }

  virtual double get_phi() const
  {
    PHOOL_VIRTUAL_WARN("get_phi()");
    return NAN;
  }

  virtual void set_tower_type(int)
  {
    PHOOL_VIRTUAL_WARN("set_tower_type()");
    return;
  }
  virtual int get_tower_type() const
  {
    PHOOL_VIRTUAL_WARN("get_tower_type()");
    return -1;
  }

  virtual double get_roty() const
  {
    PHOOL_VIRTUAL_WARN("get_roty()");
    return NAN;
  }
  virtual double get_rotz() const
  {
    PHOOL_VIRTUAL_WARN("get_rotz()");
    return NAN;
  }

  virtual void set_roty(double)
  {
    PHOOL_VIRTUAL_WARN("set_roty()");
    return;
  }
  virtual void set_rotz(double)
  {
    PHOOL_VIRTUAL_WARN("set_rotz()");
    return;
  }

 protected:
  RawTowerGeom() {}

  ClassDefOverride(RawTowerGeom, 2)
};

#endif
