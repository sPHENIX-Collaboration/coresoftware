#ifndef CALOBASE_RAWTOWERGEOMV3_H
#define CALOBASE_RAWTOWERGEOMV3_H

#include "RawTowerGeom.h"

#include "RawTowerDefs.h"

#include <cmath>
#include <iostream>

class RawTowerGeomv3 : public RawTowerGeom
{
 public:
  RawTowerGeomv3() {}
  RawTowerGeomv3(RawTowerDefs::keytype id);
  ~RawTowerGeomv3() override {}

  void identify(std::ostream& os = std::cout) const override;

  void set_id(RawTowerDefs::keytype key) override { _towerid = key; }
  RawTowerDefs::keytype get_id() const override { return _towerid; }

  int get_bineta() const override;
  int get_binphi() const override;

  int get_column() const override { return get_bineta(); }
  int get_row() const override { return get_binphi(); }

  int get_binl() const override { return RawTowerDefs::decode_index3v2(_towerid); }

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
  double get_theta() const override;

  void set_tower_type(int tt) override { _tower_type = tt; }
  int get_tower_type() const override { return _tower_type; }

 protected:
  RawTowerDefs::keytype _towerid = ~0;  // complement = 0xFFFFF... independent of integer type (32/64/... bits)

  double _center_x = NAN;
  double _center_y = NAN;
  double _center_z = NAN;

  double _size_x = NAN;
  double _size_y = NAN;
  double _size_z = NAN;

  int _tower_type = -1;

  ClassDefOverride(RawTowerGeomv3, 4)
};

#endif /* CALOBASE_RAWTOWERGEOMV3_H */
