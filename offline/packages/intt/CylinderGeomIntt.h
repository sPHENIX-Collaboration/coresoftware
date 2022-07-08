#ifndef INTT_CYLINDERGEOMINTT_H
#define INTT_CYLINDERGEOMINTT_H

#include <g4detectors/PHG4CylinderGeom.h>

#include <trackbase/ActsGeometry.h>

#include <cmath>
#include <iostream>

#include <TVector3.h>

class CylinderGeomIntt : public PHG4CylinderGeom
{
 public:
  CylinderGeomIntt();
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
  void identify(std::ostream &os = std::cout) const override;

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
  void find_segment_center(Surface surface, ActsGeometry* tGeometry, double location[]);
  void find_strip_center(Surface surface, ActsGeometry *tGeometry, const int segment_z_bin,  const int segment_phi_bin, const int strip_column, const int strip_index, double location[]);
  void find_strip_index_values(const int segment_z_bin, const double ypos, const double zpos, int &strip_y_index, int &strip_z_index) override;

  bool load_geometry();
  void find_strip_center_localcoords(const int segment_z_bin, const int strip_y_index, const int strip_z_index, double location[]);
  void find_indices_from_segment_center(int &segment_z_bin, int &segment_phi_bin, double location[]);
  TVector3 get_world_from_local_coords(Surface surface, ActsGeometry* tGeometry, TVector2 local);
  TVector3 get_world_from_local_coords(Surface surface, ActsGeometry* tGeometry, TVector3 local);
  TVector3 get_local_from_world_coords(Surface surface, ActsGeometry* tGeometry, TVector3 world);
  void find_indices_from_world_location(int &segment_z_bin, int &segment_phi_bin, double location[]);

  void find_strip_center(int, int, int, int, double*) override
  { std::cout << "find_strip_center(int, int, int, int, double[]) is deprecated" << std::endl;
  }
  void find_segment_center(const int, const int, double*) override
  { std::cout << "find_segment_center(const int, const int, double*) is deprecated" << std::endl; }

  double get_strip_phi_tilt() const
  {
    return m_OffsetRot;
  }

 protected:
  int m_Layer;
  int m_NStripsPhiCell;
  int m_NStripsZSensor[2];
  double m_StripX;
  double m_StripY;
  double m_SensorRadius;
  double m_StripXOffset;
  double m_OffsetPhi;
  double m_OffsetRot;
  double m_dPhi;

  double m_StripZ[2];
  double m_LadderZ[2];

  ClassDefOverride(CylinderGeomIntt, 1)
};

#endif
