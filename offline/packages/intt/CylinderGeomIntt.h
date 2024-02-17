#ifndef INTT_CYLINDERGEOMINTT_H
#define INTT_CYLINDERGEOMINTT_H

#include <g4detectors/PHG4CylinderGeom.h>

#include <trackbase/ActsGeometry.h>

#include <TVector2.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>

class CylinderGeomIntt : public PHG4CylinderGeom
{
 public:
  CylinderGeomIntt() = default;
  CylinderGeomIntt(
      const int layer,
      const double strip_x,
      const double strip_y,
      const double strip_z0,
      const double strip_z1,
      const int nstrips_z_sensor0,
      const int nstrips_z_sensor1,
      const int nstrips_phi_cell,
      const int nladders_layer,
      const double ladder_z0,
      const double ladder_z1,
      const double sensor_radius,
      const double strip_x_offset,
      const double offsetphi,
      const double offsetrot)
    : m_Layer(layer)
    , m_NStripsPhiCell(nstrips_phi_cell)
    , m_StripX(strip_x)
    , m_StripY(strip_y)
    , m_SensorRadius(sensor_radius)
    , m_StripXOffset(strip_x_offset)
    , m_OffsetPhi(offsetphi)
    , m_OffsetRot(offsetrot)
  {
    // Type-A
    m_StripZ[0] = strip_z0;
    m_LadderZ[0] = ladder_z0;
    m_NStripsZSensor[0] = nstrips_z_sensor0;

    // Type-B
    m_StripZ[1] = strip_z1;
    m_LadderZ[1] = ladder_z1;
    m_NStripsZSensor[1] = nstrips_z_sensor1;

    m_dPhi = 2. * M_PI / nladders_layer;
  }

  // from PHObject
  void identify(std::ostream& os = std::cout) const override;

  // overridden from base class
  double get_thickness() const override
  {
    return m_StripX;
  }

  double get_strip_y_spacing() const override
  {
    return m_StripY;
  }

  double get_strip_z_spacing() const override
  {
    return m_StripZ[0];
  }

  double get_strip_tilt() const override
  {
    return 0.;
  }

  void set_layer(const int i) override
  {
    m_Layer = i;
  }

  int get_layer() const override
  {
    return m_Layer;
  }

  double get_radius() const override
  {
    return m_SensorRadius;
  }

  // our own
  void find_segment_center(const Surface& surface, ActsGeometry* tGeometry, double location[]);
  void find_strip_center(const Surface& surface, ActsGeometry* tGeometry, const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[]);
  void find_strip_index_values(const int segment_z_bin, const double ypos, const double zpos, int& strip_y_index, int& strip_z_index) override;

  bool load_geometry() { return true; }
  void find_strip_center_localcoords(const int segment_z_bin, const int strip_y_index, const int strip_z_index, double location[]);
  void find_indices_from_segment_center(int& segment_z_bin, int& segment_phi_bin, double location[]);
  TVector3 get_world_from_local_coords(const Surface& surface, ActsGeometry* tGeometry, const TVector2& local);
  TVector3 get_world_from_local_coords(const Surface& surface, ActsGeometry* tGeometry, const TVector3& local);
  TVector3 get_local_from_world_coords(const Surface& surface, ActsGeometry* tGeometry, TVector3 world);
  void find_indices_from_world_location(int& segment_z_bin, int& segment_phi_bin, double location[]);

  void find_strip_center(int, int, int, int, double*) override
  {
    std::cout << "find_strip_center(int, int, int, int, double[]) is deprecated" << std::endl;
  }
  void find_segment_center(const int, const int, double*) override
  {
    std::cout << "find_segment_center(const int, const int, double*) is deprecated" << std::endl;
  }

  double get_strip_phi_tilt() const
  {
    return m_OffsetRot;
  }

 protected:
  int m_Layer{-1};
  int m_NStripsPhiCell{-1};
  int m_NStripsZSensor[2] {-1,-1};
  double m_StripX{std::numeric_limits<double>::quiet_NaN()};
  double m_StripY{std::numeric_limits<double>::quiet_NaN()};
  double m_SensorRadius{std::numeric_limits<double>::quiet_NaN()};
  double m_StripXOffset{std::numeric_limits<double>::quiet_NaN()};
  double m_OffsetPhi{std::numeric_limits<double>::quiet_NaN()};
  double m_OffsetRot{std::numeric_limits<double>::quiet_NaN()};
  double m_dPhi{std::numeric_limits<double>::quiet_NaN()};

  double m_StripZ[2] {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  double m_LadderZ[2]{std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};

  ClassDefOverride(CylinderGeomIntt, 1)
};

#endif
