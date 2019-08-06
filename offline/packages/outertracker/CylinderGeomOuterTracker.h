#ifndef OTRACK_CYLINDERGEOMOUTERTRACKER_H
#define OTRACK_CYLINDERGEOMOUTERTRACKER_H

#include <g4detectors/PHG4CylinderGeom.h>

#include <cmath>
#include <iostream>

class CylinderGeomOuterTracker : public PHG4CylinderGeom
{
 public:
  CylinderGeomOuterTracker();
  CylinderGeomOuterTracker(
      const int layer,
      const double inner_radius,
      const double outer_radius,
      const double length,
      const int nseg_phi,
      const int nseg_z )
     : m_Layer(layer)
    , m_NStripsPhi(nseg_phi)
    , m_NStripsZ(nseg_z)
    , m_InnerRadius(inner_radius)
    , m_OuterRadius(outer_radius)
    , m_Length(length)
  {
    m_dPhi = 2.0 * M_PI / (double) m_NStripsPhi;
    m_dZ = m_Length / (double) m_NStripsZ;

    return;
  }

  void identify(std::ostream &os = std::cout) const;
  void set_layer(const int i)
  {
    m_Layer = i;
  }

  int get_layer() const
  {
    return m_Layer;
  }

  double get_inner_radius() const
  {
    return m_InnerRadius;
  }

  double get_outer_radius() const
  {
    return m_OuterRadius;
  }

  double get_average_radius() const
  {
    return (m_OuterRadius + m_InnerRadius) / 2.0;
  }

  bool load_geometry();
  void find_pixel_center(const int pixel_z_bin, const int pixel_phi_bin, double &phi, double &z);
  void find_pixel_index_values(const double xpos, const double ypos, const double zpos, int &pixel_phi_index, int &pixel_z_index);

  double get_thickness() const
  {
    return m_OuterRadius - m_InnerRadius;
  }

  double get_pixel_phi_spacing() const
  {
    return m_dPhi;
  }

  double get_pixel_z_spacing() const
  {
    return m_dZ;
  }

  int get_num_pixels_z() const
  {
    return m_NStripsZ;
  }

  int get_num_pixels_phi() const
  {
    return m_NStripsPhi;
  }

 protected:
  int m_Layer;
  int m_NStripsPhi;
  int m_NStripsZ;
  double m_InnerRadius;
  double m_OuterRadius;
  double m_Length;
  double m_dZ;
  double m_dPhi;

  ClassDef(CylinderGeomOuterTracker, 1)
};

#endif
