#include "CylinderGeomOuterTracker.h"

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/Rotation.h>

#include <algorithm>
#include <cmath>

using namespace std;

CylinderGeomOuterTracker::CylinderGeomOuterTracker()
  : m_Layer(-1)
  , m_NStripsPhi(-1)
  , m_NStripsZ(-1)
  , m_InnerRadius(NAN)
  , m_OuterRadius(NAN)
  , m_Length(NAN)
  , m_dZ(NAN)
  , m_dPhi(NAN)
{
  return;
}

void CylinderGeomOuterTracker::identify(std::ostream &os) const
{
  os << "CylinderGeomOuterTracker Object" << endl;
  os << "layer: " << get_layer() << endl;
  os << "InnerRadius: " << get_inner_radius() << endl;
  os << "OuterRadius: " << get_outer_radius() << endl;
  os << " Phi segments " << m_NStripsPhi << endl;
  os << " Z segments " << m_NStripsZ << endl;
}

bool CylinderGeomOuterTracker::load_geometry()
{
  return true;
}

void CylinderGeomOuterTracker::find_pixel_center(const int pixel_z_bin, const int pixel_phi_bin, double &phi, double &z)
{
  //double pixel_radius = (m_InnerRadius + m_OuterRadius) / 2.0;
  // Ladder
  phi = m_dPhi *  pixel_phi_bin;
  z = m_dZ * pixel_z_bin;
  /*
  location[0] = pixel_radius * cos(phi);
  location[1] = pixel_radius * sin(phi);
  location[2] = m_dZ * pixel_z_bin;
  */
}

void CylinderGeomOuterTracker::find_pixel_index_values(const double xin, const double yin, const double zin, int &pixel_y_index, int &pixel_z_index)
{
  // expect cm

  // get the strip z index
  pixel_z_index = (int) (zin / m_dZ);

  //double sensor_radius = (m_InnerRadius + m_OuterRadius) / 2.0;
  double phi_in = atan(yin/xin);
  // get phi_in in the correct sector
  cout << " Fix quadrant: phi_in " << phi_in << " yin " << yin << " xin " << xin << endl;

  // get the strip y index
  pixel_y_index = (int) (phi_in / m_dPhi);

}

