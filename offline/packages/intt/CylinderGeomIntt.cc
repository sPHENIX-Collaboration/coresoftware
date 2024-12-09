#include "CylinderGeomIntt.h"

#include <algorithm>
#include <cmath>
#include <memory>  // for __shared_ptr_access
#include <utility>

void CylinderGeomIntt::identify(std::ostream& os) const
{
  os << "CylinderGeomIntt Object" << std::endl;
  os << "layer: " << get_layer() << std::endl;
  os << "Radius: " << get_radius() << std::endl;
}

void CylinderGeomIntt::find_indices_from_world_location(int& segment_z_bin, int& segment_phi_bin, double location[])
{
  double signz = (location[2] > 0) ? 1. : -1;
  double phi = atan2(location[1], location[0]);
  double tolerance_phi = 0.05;
  double tolerance_z = 0.5;

  if (fabs(phi - m_OffsetPhi) > tolerance_phi && phi < 0)
  {
    phi += 2.0 * M_PI;
  }
  double segment_phi_bin_tmp = (phi - m_OffsetPhi) / m_dPhi;
  segment_phi_bin = round(segment_phi_bin_tmp);

  double z_tmp = location[2] / signz;

  // decide if this is a type A (0) or type B (1) sensor
  int itype;
  // if (fabs((z_tmp / m_LadderZ[0])) < 1.0)
  if (fabs((1.0 - z_tmp / m_LadderZ[0])) < tolerance_z)
  {
    itype = 0;
  }
  else
  {
    itype = 1;
  }

  if (signz < 0)
  {
    segment_z_bin = itype;  // 0 = itype 0 +z,  1 = itype 1 +z,  2 = itupe 0 -z, 3 = itype 1 -z
  }
  else
  {
    segment_z_bin = itype + 2;
  }
}

void CylinderGeomIntt::find_indices_from_segment_center(int& segment_z_bin, int& segment_phi_bin, double location[])
{
  double signz = (location[2] > 0) ? 1. : -1;
  double phi = atan2(location[1], location[0]);
  // std::cout << "phi before 2pi shift=" << phi << " offset " << m_OffsetPhi << " fabs(phi - m_OffsetPhi)=" << fabs(phi - m_OffsetPhi) << std::endl;
  double tolerance_phi = 0.05;
  double tolerance_z = 0.5;
  if (fabs(phi - m_OffsetPhi) > tolerance_phi && phi < 0)
  {
    phi += 2.0 * M_PI;
  }
  double segment_phi_bin_tmp = (phi - m_OffsetPhi) / m_dPhi;
  segment_phi_bin = lround(segment_phi_bin_tmp);

  // std::cout << "     phi " <<phi << " segment_phi_bin_tmp " <<  segment_phi_bin_tmp << " segment_phi_bin " << segment_phi_bin << " location " << location[0] << "  " << location[1] << "  " << location[2] << std::endl;

  double z_tmp = location[2] / signz;

  // decide if this is a type A (0) or type B (1) sensor
  int itype;
  if (fabs((1.0 - z_tmp / m_LadderZ[0])) < tolerance_z)
  {
    itype = 0;
  }
  else
  {
    itype = 1;
  }

  if (signz < 0)
  {
    segment_z_bin = itype;  // 0 = itype 0 +z,  1 = itype 1 +z,  2 = itupe 1 -z, 3 = itype 1 -z
  }
  else
  {
    segment_z_bin = itype + 2;
  }

  // std::cout << " world coords: " <<  location[0] << " " << location[1] << " " << location[2] <<  " signz " << signz << " itype " << itype << " z_tmp " << z_tmp <<  " m_LadderZ " << m_LadderZ[itype] << std::endl;
  // std::cout << "radius " << m_SensorRadius << " offsetphi " << m_OffsetPhi << " rad  dphi_ " << m_dPhi << " rad  segment_phi_bin " << segment_phi_bin << " phi " << phi  << std::endl;
}

void CylinderGeomIntt::find_strip_index_values(const int segment_z_bin, const double yin, const double zin, int& strip_y_index, int& strip_z_index)
{
  // Given the location in y and z in sensor local coordinates, find the strip y and z index values

  // find the sensor type (inner or outer) from the segment_z_bin (location of sensor on ladder)
  const int itype = segment_z_bin % 2;
  if (itype != 0 && itype != 1)
  {
    std::cout << "Problem: itype = " << itype << std::endl;
    return;
  }

  // expect cm
  double zpos = zin;
  double ypos = yin;

  const double strip_z = m_StripZ[itype];
  const int nstrips_z_sensor = m_NStripsZSensor[itype];
  const int nstrips_y_sensor = m_NStripsPhiCell;

  // get the strip z index
  double zup = (double) nstrips_z_sensor * strip_z / 2.0 + zpos;
  strip_z_index = floor(zup / strip_z);

  // get the strip y index
  double yup = (double) nstrips_y_sensor * m_StripY / 2.0 + ypos;
  strip_y_index = floor(yup / m_StripY);

  /*
  std::cout << "segment_z_bin " << segment_z_bin << " ypos " << ypos << " zpos " << zpos << " zup " << zup << " yup " << yup << std::endl;
  std::cout << "      -- itype " << itype << " strip_y " << m_StripY << " strip_z " << strip_z << " nstrips_z_sensor " << nstrips_z_sensor
       << " nstrips_y_sensor " << nstrips_y_sensor << std::endl;
  std::cout << "      --  strip_z_index " << strip_z_index << " strip_y_index " << strip_y_index << std::endl;
  */
}

// this name is a really bad idea whcih ticks off clang-tidy (justifiably) but it seems intentional
// NOLINTNEXTLINE(bugprone-virtual-near-miss)
void CylinderGeomIntt::find_strip_center_localcoords(const int segment_z_bin, const int strip_y_index, const int strip_z_index, double location[])
{
  // find the sensor type (inner or outer) from the segment_z_bin (location of sensor on ladder)
  const int itype = segment_z_bin % 2;
  if (itype != 0 && itype != 1)
  {
    std::cout << "Problem: itype = " << itype << std::endl;
    return;
  }

  const double strip_z = m_StripZ[itype];
  const int nstrips_z_sensor = m_NStripsZSensor[itype];
  const int nstrips_y_sensor = m_NStripsPhiCell;

  // center of strip in y
  double ypos = (double) strip_y_index * m_StripY + m_StripY / 2.0 - (double) nstrips_y_sensor * m_StripY / 2.0;

  // center of strip in z
  double zpos = (double) strip_z_index * strip_z + strip_z / 2.0 - (double) nstrips_z_sensor * strip_z / 2.0;

  location[0] = 0.0;
  location[1] = ypos;
  location[2] = zpos;
}
