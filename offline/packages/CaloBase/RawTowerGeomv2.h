#ifndef CALOBASE_RAWTOWERGEOMV2_H
#define CALOBASE_RAWTOWERGEOMV2_H

#include "RawTowerGeom.h"

#include "RawTowerDefs.h"

#include <iostream>

class RawTowerGeomv2 : public RawTowerGeom
{
 public:
  RawTowerGeomv2() {}
  RawTowerGeomv2(RawTowerDefs::keytype id);
  ~RawTowerGeomv2() override {}

  void identify(std::ostream& os = std::cout) const override;

  void set_id(RawTowerDefs::keytype key) override { _towerid = key; }
  RawTowerDefs::keytype get_id() const override { return _towerid; }

  int get_bineta() const override { return RawTowerDefs::decode_index1(_towerid); }
  int get_binphi() const override { return RawTowerDefs::decode_index2(_towerid); }
  int get_column() const override { return RawTowerDefs::decode_index1(_towerid); }
  int get_row() const override { return RawTowerDefs::decode_index2(_towerid); }

  void set_center_x(double x) override
  {
    _center_x = x;
    return;
  }
  void set_center_y(double y) override
  {
    _center_y = y;
    return;
  }
  void set_center_z(double z) override
  {
    _center_z = z;
    return;
  }

  void set_size_x(double dx) override
  {
    _size_x = dx;
    return;
  }
  void set_size_y(double dy) override
  {
    _size_y = dy;
    return;
  }
  void set_size_z(double dz) override
  {
    _size_z = dz;
    return;
  }

  double get_center_x() const override { return _center_x; }
  double get_center_y() const override { return _center_y; }
  double get_center_z() const override { return _center_z; }

  double get_size_x() const override { return _size_x; }
  double get_size_y() const override { return _size_y; }
  double get_size_z() const override { return _size_z; }
  double get_volume() const override { return (_size_x * _size_y * _size_z); }

  double get_center_radius() const override;
  double get_eta() const override;
  double get_phi() const override;

 protected:
  RawTowerDefs::keytype _towerid = ~0U;

  double _center_x = 0.;
  double _center_y = 0.;
  double _center_z = 0.;

  double _size_x = 0.;
  double _size_y = 0.;
  double _size_z = 0.;

  ClassDefOverride(RawTowerGeomv2, 3)
};

#endif
