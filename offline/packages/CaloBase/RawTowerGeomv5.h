#ifndef CALOBASE_RAWTOWERGEOMV5_H
#define CALOBASE_RAWTOWERGEOMV5_H

#include "RawTowerGeom.h"

#include "RawTowerDefs.h"

#include <iostream>
#include <limits>

class RawTowerGeomv5 : public RawTowerGeom
{
 public:

  RawTowerGeomv5() {}
  RawTowerGeomv5(RawTowerDefs::keytype id);
  RawTowerGeomv5(const RawTowerGeom& geom0);
  ~RawTowerGeomv5() override {}

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
  void set_rotx(double rotx) override
  {
    _rotx = rotx;
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

  void set_vertices(const std::vector<double>&) override;

  double get_center_x() const override { return _center_x; }
  double get_center_y() const override { return _center_y; }
  double get_center_z() const override { return _center_z; }
  double get_center_int_x() const override { return _center_int_x; }
  double get_center_int_y() const override { return _center_int_y; }
  double get_center_int_z() const override { return _center_int_z; }
  double get_center_ext_x() const override { return _center_ext_x; }
  double get_center_ext_y() const override { return _center_ext_y; }
  double get_center_ext_z() const override { return _center_ext_z; }
  double get_center_low_eta_x() const override { return _center_low_eta_x; }
  double get_center_low_eta_y() const override { return _center_low_eta_y; }
  double get_center_low_eta_z() const override { return _center_low_eta_z; }
  double get_center_high_eta_x() const override { return _center_high_eta_x; }
  double get_center_high_eta_y() const override { return _center_high_eta_y; }
  double get_center_high_eta_z() const override { return _center_high_eta_z; }
  double get_center_low_phi_x() const override { return _center_low_phi_x; }
  double get_center_low_phi_y() const override { return _center_low_phi_y; }
  double get_center_low_phi_z() const override { return _center_low_phi_z; }
  double get_center_high_phi_x() const override { return _center_high_phi_x; }
  double get_center_high_phi_y() const override { return _center_high_phi_y; }
  double get_center_high_phi_z() const override { return _center_high_phi_z; }
  double get_rotx() const override { return _rotx; }
  double get_roty() const override { return _roty; }
  double get_rotz() const override { return _rotz; }
  double get_vertex_x(int i) const override { return _vertices_x[i]; }
  double get_vertex_y(int i) const override { return _vertices_y[i]; }
  double get_vertex_z(int i) const override { return _vertices_z[i]; }
  double get_center_radius() const override;
  double get_eta() const override;
  double get_phi() const override;
  double get_theta() const override;

 protected:
  
  RawTowerDefs::keytype _towerid{std::numeric_limits<RawTowerDefs::keytype>::max()};

  static constexpr int _nVtx = 8;
  double _center_x{std::numeric_limits<double>:: quiet_NaN()};
  double _center_y{std::numeric_limits<double>:: quiet_NaN()};
  double _center_z{std::numeric_limits<double>:: quiet_NaN()};
  double _center_int_x{std::numeric_limits<double>:: quiet_NaN()};
  double _center_int_y{std::numeric_limits<double>:: quiet_NaN()};
  double _center_int_z{std::numeric_limits<double>:: quiet_NaN()};
  double _center_ext_x{std::numeric_limits<double>:: quiet_NaN()};
  double _center_ext_y{std::numeric_limits<double>:: quiet_NaN()};
  double _center_ext_z{std::numeric_limits<double>:: quiet_NaN()};
  double _center_low_eta_x{std::numeric_limits<double>:: quiet_NaN()};
  double _center_low_eta_y{std::numeric_limits<double>:: quiet_NaN()};
  double _center_low_eta_z{std::numeric_limits<double>:: quiet_NaN()};
  double _center_high_eta_x{std::numeric_limits<double>:: quiet_NaN()};
  double _center_high_eta_y{std::numeric_limits<double>:: quiet_NaN()};
  double _center_high_eta_z{std::numeric_limits<double>:: quiet_NaN()};
  double _center_low_phi_x{std::numeric_limits<double>:: quiet_NaN()};
  double _center_low_phi_y{std::numeric_limits<double>:: quiet_NaN()};
  double _center_low_phi_z{std::numeric_limits<double>:: quiet_NaN()};
  double _center_high_phi_x{std::numeric_limits<double>:: quiet_NaN()};
  double _center_high_phi_y{std::numeric_limits<double>:: quiet_NaN()};
  double _center_high_phi_z{std::numeric_limits<double>:: quiet_NaN()};

  double _vertices_x[_nVtx] = {std::numeric_limits<double>:: quiet_NaN()};
  double _vertices_y[_nVtx] = {std::numeric_limits<double>:: quiet_NaN()};
  double _vertices_z[_nVtx] = {std::numeric_limits<double>:: quiet_NaN()};

  double _rotx{std::numeric_limits<double>:: quiet_NaN()};
  double _roty{std::numeric_limits<double>:: quiet_NaN()};
  double _rotz{std::numeric_limits<double>:: quiet_NaN()};
  
  ClassDefOverride(RawTowerGeomv5, 4)
};

#endif
