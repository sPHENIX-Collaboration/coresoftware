#ifndef CALOBASE_RAWTOWERGEOMV1_H
#define CALOBASE_RAWTOWERGEOMV1_H

#include "RawTowerGeom.h"

#include "RawTowerDefs.h"

#include <cmath>
#include <iostream>

class RawTowerGeomv1 : public RawTowerGeom
{
 public:
  RawTowerGeomv1() {}
  RawTowerGeomv1(RawTowerDefs::keytype id);
  ~RawTowerGeomv1() override {}

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

  double get_center_x() const override { return _center_x; }
  double get_center_y() const override { return _center_y; }
  double get_center_z() const override { return _center_z; }

  double get_center_radius() const override;
  double get_eta() const override;
  double get_phi() const override;
  double get_theta() const override;

 protected:
  RawTowerDefs::keytype _towerid = ~0;  // 0xFFFFFF.. independant of type

  double _center_x = NAN;
  double _center_y = NAN;
  double _center_z = NAN;
  ClassDefOverride(RawTowerGeomv1, 4)
};

#endif
