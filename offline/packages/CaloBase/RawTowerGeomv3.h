#ifndef CALOBASE_RAWTOWERGEOMV3_H
#define CALOBASE_RAWTOWERGEOMV3_H

#include "RawTowerGeom.h"

class RawTowerGeomv3 : public RawTowerGeom
{
 public:
  RawTowerGeomv3();
  RawTowerGeomv3(RawTowerDefs::keytype id);
  virtual ~RawTowerGeomv3() {}

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

  void set_size_x(double dx)
  {
    _size_x = dx;
    return;
  }
  void set_size_y(double dy)
  {
    _size_y = dy;
    return;
  }
  void set_size_z(double dz)
  {
    _size_z = dz;
    return;
  }

  double get_center_x() const { return _center_x; }
  double get_center_y() const { return _center_y; }
  double get_center_z() const { return _center_z; }

  double get_size_x() const { return _size_x; }
  double get_size_y() const { return _size_y; }
  double get_size_z() const { return _size_z; }
  double get_volume() const { return (_size_x * _size_y * _size_z); }

  double get_center_radius() const;
  double get_eta() const;
  double get_phi() const;

  void set_tower_type(int tt)
  {
    _tower_type = tt;
    return;
  }
  int get_tower_type() const { return _tower_type; }

 protected:
  RawTowerDefs::keytype _towerid;

  double _center_x;
  double _center_y;
  double _center_z;

  double _size_x;
  double _size_y;
  double _size_z;

  int _tower_type;

  ClassDef(RawTowerGeomv3, 4)
};

#endif /* CALOBASE_RAWTOWERGEOMV3_H */
