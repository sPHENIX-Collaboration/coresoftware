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
  virtual ~RawTowerGeomv1() {}

  void identify(std::ostream& os = std::cout) const;

  void set_id(RawTowerDefs::keytype key) { _towerid = key; }
  RawTowerDefs::keytype get_id() const { return _towerid; }

  int get_bineta() const { return RawTowerDefs::decode_index1(_towerid); }
  int get_binphi() const { return RawTowerDefs::decode_index2(_towerid); }
  int get_column() const { return RawTowerDefs::decode_index1(_towerid); }
  int get_row() const { return RawTowerDefs::decode_index2(_towerid); }

  void set_center_x(double x)
  {
    _center_x = x;
    return;
  }
  void set_center_y(double y)
  {
    _center_y = y;
    return;
  }
  void set_center_z(double z)
  {
    _center_z = z;
    return;
  }

  double get_center_x() const { return _center_x; }
  double get_center_y() const { return _center_y; }
  double get_center_z() const { return _center_z; }

  double get_center_radius() const;
  double get_eta() const;
  double get_phi() const;
  double get_theta() const;

 protected:
  RawTowerDefs::keytype _towerid = ~0;  // 0xFFFFFF.. independant of type

  double _center_x = NAN;
  double _center_y = NAN;
  double _center_z = NAN;
  ClassDef(RawTowerGeomv1, 4)
};

#endif
