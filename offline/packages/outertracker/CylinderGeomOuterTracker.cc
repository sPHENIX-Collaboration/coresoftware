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

void CylinderGeomOuterTracker::find_pixel_center(const int pixel_phi_bin, const int pixel_z_bin, double &phi, double &z)
{
  phi = m_dPhi / 2.0 + m_dPhi *  pixel_phi_bin - M_PI;
  z = m_dZ / 2.0 + m_dZ * pixel_z_bin - m_Length/2.0;

  //cout << " pixel_phi_bin " << pixel_phi_bin << " m_dPhi " << m_dPhi << " phi " << phi << endl;
  //cout << " pixel_z_bin " << pixel_z_bin << " m_dZ " << m_dZ << " z " << z << endl;
}

void CylinderGeomOuterTracker::find_pixel_index_values(const double xin, const double yin, const double zin, int &pixel_y_index, int &pixel_z_index)
{
  // expect cm

  // get the strip z index - z is from -m_Length/2 to + m_Length/2
  pixel_z_index = (int) ( (zin + m_Length/2.0) / m_dZ);

  double phi_in = atan2(yin, xin);

  // get the strip y index -   phi is from -pi to +pi
  pixel_y_index = (int) ( (phi_in + M_PI ) / m_dPhi);

  //cout << " phi_in " << phi_in << " yin " << yin << " xin " << xin << " zin " << zin << endl;
  //cout << " pixel y index " << pixel_y_index << " pixel_z_index " << pixel_z_index << endl;

}

