#ifndef CALOBASE_RAWTOWERGEOMV5_H
#define CALOBASE_RAWTOWERGEOMV5_H

#include "RawTowerGeom.h"

#include "RawTowerDefs.h"

#include <iostream>
#include <limits>

class RawTowerGeomv5 : public RawTowerGeom
{
 public:
  RawTowerGeomv5() = default;
  explicit RawTowerGeomv5(RawTowerDefs::keytype id);
  explicit RawTowerGeomv5(const RawTowerGeom& geom0);
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
    _center.x = x;
    return;
  }
  void set_center_y(double y) override
  {
    _center.y = y;
    return;
  }
  void set_center_z(double z) override
  {
    _center.z = z;
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

  void set_vertices(const std::vector<double>& /*vertices*/) override;

  double get_center_x() const override { return _center.x; }
  double get_center_y() const override { return _center.y; }
  double get_center_z() const override { return _center.z; }
  double get_center_int_x() const override { return _center_int.x; }
  double get_center_int_y() const override { return _center_int.y; }
  double get_center_int_z() const override { return _center_int.z; }
  double get_center_ext_x() const override { return _center_ext.x; }
  double get_center_ext_y() const override { return _center_ext.y; }
  double get_center_ext_z() const override { return _center_ext.z; }
  double get_center_low_eta_x() const override { return _center_low_eta.x; }
  double get_center_low_eta_y() const override { return _center_low_eta.y; }
  double get_center_low_eta_z() const override { return _center_low_eta.z; }
  double get_center_high_eta_x() const override { return _center_high_eta.x; }
  double get_center_high_eta_y() const override { return _center_high_eta.y; }
  double get_center_high_eta_z() const override { return _center_high_eta.z; }
  double get_center_low_phi_x() const override { return _center_low_phi.x; }
  double get_center_low_phi_y() const override { return _center_low_phi.y; }
  double get_center_low_phi_z() const override { return _center_low_phi.z; }
  double get_center_high_phi_x() const override { return _center_high_phi.x; }
  double get_center_high_phi_y() const override { return _center_high_phi.y; }
  double get_center_high_phi_z() const override { return _center_high_phi.z; }
  double get_rotx() const override { return _rotx; }
  double get_roty() const override { return _roty; }
  double get_rotz() const override { return _rotz; }
  double get_vertex_x(int i) const override { return _vertices[i].x; }
  double get_vertex_y(int i) const override { return _vertices[i].y; }
  double get_vertex_z(int i) const override { return _vertices[i].z; }
  double get_center_radius() const override;
  double get_eta() const override;
  double get_phi() const override;
  double get_theta() const override;

  // This structure is only used by this class
  struct point_coordinates
  {
    double x;
    double y;
    double z;
  };

 protected:
  using vertex_t = std::pair<double, double>;
  using vertices_t = std::vector<vertex_t>;

  RawTowerDefs::keytype _towerid{std::numeric_limits<RawTowerDefs::keytype>::max()};

  static constexpr int _nVtx = 8;
  point_coordinates _center{std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::quiet_NaN()};
  point_coordinates _center_int{std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN()};
  point_coordinates _center_ext{std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN()};
  point_coordinates _center_low_eta{std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::quiet_NaN()};
  point_coordinates _center_high_eta{std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN()};
  point_coordinates _center_low_phi{std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::quiet_NaN()};
  point_coordinates _center_high_phi{std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN()};

  std::vector<point_coordinates> _vertices;
  double _rotx{std::numeric_limits<double>::quiet_NaN()};
  double _roty{std::numeric_limits<double>::quiet_NaN()};
  double _rotz{std::numeric_limits<double>::quiet_NaN()};

  ClassDefOverride(RawTowerGeomv5, 4)
};

#endif
