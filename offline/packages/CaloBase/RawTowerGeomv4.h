#ifndef CALOBASE_RAWTOWERGEOMV4_H
#define CALOBASE_RAWTOWERGEOMV4_H

#include "RawTowerGeom.h"

#include "RawTowerDefs.h"

#include <TRotation.h>
#include <TVector3.h>
#include <cmath>
#include <iostream>

class RawTowerGeomv4 : public RawTowerGeom
{
 public:
  RawTowerGeomv4() {}
  RawTowerGeomv4(RawTowerDefs::keytype id);
  virtual ~RawTowerGeomv4() {}

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
  void set_roty(double roty) override
  {
    _roty = roty;
    return;
  }
  void set_rotz(double rotz) override
  {
    _rotz = rotz;
    return;
  }

  double get_center_x() const override { return _center_x; }
  double get_center_y() const override { return _center_y; }
  double get_center_z() const override { return _center_z; }
  double get_roty() const override { return _roty; }
  double get_rotz() const override { return _rotz; }

  TVector3 get_final_position() const;
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
  double _roty = NAN;
  double _rotz = NAN;

  int _tower_type = -1;

  ClassDefOverride(RawTowerGeomv4, 5)
};

#endif /* CALOBASE_RawTowerGeomv4_H */
