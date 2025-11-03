#ifndef CALOBASE_RAWTOWERGEOMV5_H
#define CALOBASE_RAWTOWERGEOMV5_H

#include "RawTowerGeom.h"

#include "RawTowerDefs.h"

#include <array>
#include <iostream>
#include <limits>

class RawTowerGeomv5 : public RawTowerGeom
{
 public:
  RawTowerGeomv5() = default;
  RawTowerGeomv5(RawTowerDefs::keytype id);
  RawTowerGeomv5(const RawTowerGeom& geom0);
  ~RawTowerGeomv5() override = default;

  void identify(std::ostream& os = std::cout) const override;

  void set_id(RawTowerDefs::keytype key) override { _towerid = key; }
  RawTowerDefs::keytype get_id() const override { return _towerid; }

  int get_bineta() const override { return RawTowerDefs::decode_index1(_towerid); }
  int get_binphi() const override { return RawTowerDefs::decode_index2(_towerid); }
  int get_column() const override { return RawTowerDefs::decode_index1(_towerid); }
  int get_row() const override { return RawTowerDefs::decode_index2(_towerid); }

  void set_center_x(double x) override
  {
    _center[0] = x;
    return;
  }
  void set_center_y(double y) override
  {
    _center[1] = y;
    return;
  }
  void set_center_z(double z) override
  {
    _center[2] = z;
    return;
  }
  void set_rotx(double rotx) override
  {
    _rot[0] = rotx;
    return;
  }
  void set_roty(double roty) override
  {
    _rot[1] = roty;
    return;
  }
  void set_rotz(double rotz) override
  {
    _rot[2] = rotz;
    return;
  }

  void set_vertices(const std::vector<double>&) override;

  double get_center_x() const override { return _center[0]; }
  double get_center_y() const override { return _center[1]; }
  double get_center_z() const override { return _center[2]; }
  double get_center_int_x() const override { return _center_int[0]; }
  double get_center_int_y() const override { return _center_int[1]; }
  double get_center_int_z() const override { return _center_int[2]; }
  double get_center_ext_x() const override { return _center_ext[0]; }
  double get_center_ext_y() const override { return _center_ext[1]; }
  double get_center_ext_z() const override { return _center_ext[2]; }
  double get_center_low_eta_x() const override { return _center_low_eta[0]; }
  double get_center_low_eta_y() const override { return _center_low_eta[1]; }
  double get_center_low_eta_z() const override { return _center_low_eta[2]; }
  double get_center_high_eta_x() const override { return _center_high_eta[0]; }
  double get_center_high_eta_y() const override { return _center_high_eta[1]; }
  double get_center_high_eta_z() const override { return _center_high_eta[2]; }
  double get_center_low_phi_x() const override { return _center_low_phi[0]; }
  double get_center_low_phi_y() const override { return _center_low_phi[1]; }
  double get_center_low_phi_z() const override { return _center_low_phi[2]; }
  double get_center_high_phi_x() const override { return _center_high_phi[0]; }
  double get_center_high_phi_y() const override { return _center_high_phi[1]; }
  double get_center_high_phi_z() const override { return _center_high_phi[2]; }
  double get_rotx() const override { return _rot[0]; }
  double get_roty() const override { return _rot[1]; }
  double get_rotz() const override { return _rot[2]; }
  double get_vertex_x(int i) const override;
  double get_vertex_y(int i) const override;
  double get_vertex_z(int i) const override;
  double get_center_radius() const override;
  double get_eta() const override;
  double get_phi() const override;
  double get_theta() const override;

 protected:
  RawTowerDefs::keytype _towerid{std::numeric_limits<RawTowerDefs::keytype>::max()};

  static constexpr int _nVtx = 8;
  static constexpr int _nDim = 3;

  std::array<double, _nDim> _center{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  std::array<double, _nDim> _center_int{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  std::array<double, _nDim> _center_ext{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  std::array<double, _nDim> _center_low_eta{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  std::array<double, _nDim> _center_high_eta{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  std::array<double, _nDim> _center_low_phi{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};
  std::array<double, _nDim> _center_high_phi{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};

  std::array<double, _nVtx * _nDim> _vertices{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),  // vertex 1 (new line for readability)

      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),  // vertex 2

      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),  // etc.

      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),

      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),

      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),

      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),

      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};

  std::array<double, _nDim> _rot{
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN()};

  ClassDefOverride(RawTowerGeomv5, 1)
};

#endif
