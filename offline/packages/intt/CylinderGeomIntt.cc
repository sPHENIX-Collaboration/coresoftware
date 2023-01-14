#include "CylinderGeomIntt.h"

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/Rotation.h>

#include <Acts/Definitions/Units.hpp>
#include <algorithm>
#include <cmath>

using namespace std;

CylinderGeomIntt::CylinderGeomIntt()
  : m_Layer(-1)
  , m_NStripsPhiCell(-1)
  , m_StripX(NAN)
  , m_StripY(NAN)
  , m_SensorRadius(NAN)
  , m_StripXOffset(NAN)
  , m_OffsetPhi(NAN)
  , m_OffsetRot(NAN)
  , m_dPhi(NAN)
{
  fill_n(m_StripZ, sizeof(m_StripZ) / sizeof(double), NAN);
  fill_n(m_LadderZ, sizeof(m_LadderZ) / sizeof(double), NAN);
  fill_n(m_NStripsZSensor, sizeof(m_NStripsZSensor), -1);
  return;
}

void CylinderGeomIntt::identify(std::ostream &os) const
{
  os << "CylinderGeomIntt Object" << endl;
  os << "layer: " << get_layer() << endl;
  os << "Radius: " << get_radius() << endl;
}

bool CylinderGeomIntt::load_geometry()
{
  return true;
}

TVector3 CylinderGeomIntt::get_world_from_local_coords(Surface surface, ActsGeometry* tGeometry, TVector3 local)
{
  Acts::Vector3 loc(local.x(), local.y(), local.z());
  loc *= Acts::UnitConstants::cm;
  
  Acts::Vector3 glob = surface->transform(tGeometry->geometry().getGeoContext()) * loc;
  glob /= Acts::UnitConstants::cm;
  return TVector3(glob(0), glob(1), glob(2));
}
TVector3 CylinderGeomIntt::get_world_from_local_coords(Surface surface, ActsGeometry* tGeometry, TVector2 local)
{
  Acts::Vector2 actslocal;
  actslocal(0) = local.X();
  actslocal(1) = local.Y();
  actslocal *= Acts::UnitConstants::cm;
  
  /// Acts requires a dummy vector to be passed in the arg list
  auto global = surface->localToGlobal(tGeometry->geometry().getGeoContext(),
				       actslocal,
				       Acts::Vector3(1,1,1));

  global /= Acts::UnitConstants::cm;
  
  TVector3 ret;
  ret[0] = global(0);
  ret[1] = global(1);
  ret[2] = global(2);
  
  return ret;

}
TVector3 CylinderGeomIntt::get_local_from_world_coords(Surface surface, ActsGeometry* tGeometry, TVector3 world)
{
  Acts::Vector3 global;
  global(0) = world[0];
  global(1) = world[1];
  global(2) = world[2];
  global *= Acts::UnitConstants::cm;

  Acts::Vector3 local = surface->transform(tGeometry->geometry().getGeoContext()).inverse() * global;

  local /= Acts::UnitConstants::cm;

  /// The acts transform is offset by one element
  return TVector3(local(2), local(0), local(1));
}

void CylinderGeomIntt::find_segment_center(Surface surface, ActsGeometry* tGeometry, double location[])
{
  TVector2 local(0.0,0.0);

  TVector3 global = get_world_from_local_coords(surface, tGeometry, local);
  location[0] = global.X();
  location[1] = global.Y();
  location[2] = global.Z();
  return;
}

void CylinderGeomIntt::find_indices_from_world_location(int &segment_z_bin, int &segment_phi_bin, double location[])
{
  double signz = (location[2] > 0)? 1.  : -1;
  double phi = atan2(location[1], location[0]);
  if(phi < 0) phi += 2.0*M_PI;
  double segment_phi_bin_tmp = (phi - m_OffsetPhi)/m_dPhi;
  segment_phi_bin = round(segment_phi_bin_tmp);

  double z_tmp = location[2]  / signz;

  // decide if this is a type A (0) or type B (1) sensor
  int itype;
  if( fabs((z_tmp / m_LadderZ[0])) < 1.0) 
    itype = 0;
  else
    itype = 1;

  if(signz <0)
    segment_z_bin = itype;    // 0 = itype 0 +z,  1 = itype 1 +z,  2 = itupe 0 -z, 3 = itype 1 -z
  else
    segment_z_bin = itype + 2;
}

void CylinderGeomIntt::find_indices_from_segment_center(int &segment_z_bin, int &segment_phi_bin, double location[])
{
  double signz = (location[2] > 0)? 1.  : -1;
  double phi = atan2(location[1], location[0]);
  if(phi < 0) phi += 2.0*M_PI;
  double segment_phi_bin_tmp = (phi - m_OffsetPhi)/m_dPhi;
  segment_phi_bin = lround(segment_phi_bin_tmp);

  double z_tmp = location[2]  / signz;

  // decide if this is a type A (0) or type B (1) sensor
  int itype;
  if( fabs((1.0 - z_tmp / m_LadderZ[0])) < 0.01) 
    itype = 0;
  else
    itype = 1;

  if(signz <0)
    segment_z_bin = itype;    // 0 = itype 0 +z,  1 = itype 1 +z,  2 = itupe 1 -z, 3 = itype 1 -z
  else
    segment_z_bin = itype + 2;
  
  //cout << " world coords: " <<  location[0] << " " << location[1] << " " << location[2] <<  " signz " << signz << " itype " << itype << " z_tmp " << z_tmp <<  " m_LadderZ " << m_LadderZ[itype] << endl;
  //cout << "radius " << m_SensorRadius << " offsetphi " << m_OffsetPhi << " rad  dphi_ " << m_dPhi << " rad  segment_phi_bin " << segment_phi_bin << " phi " << phi  << endl;
}

void CylinderGeomIntt::find_strip_center(Surface surface, ActsGeometry *tGeometry, const int segment_z_bin, const int segment_phi_bin, const int strip_column, const int strip_index, double location[])
{
  // Ladder
  find_segment_center(surface, tGeometry, location);
  CLHEP::Hep3Vector ladder(location[0], location[1], location[2]);

  // Strip
  const int itype = segment_z_bin % 2;
  const double strip_z = m_StripZ[itype];
  const int nstrips_z_sensor = m_NStripsZSensor[itype];

  const double strip_localpos_z = strip_z * (strip_column % nstrips_z_sensor) - strip_z / 2. * nstrips_z_sensor + strip_z / 2.;
  // distance from bottom of sensor = m_StripY*strip_index +m_StripY/2.0, then subtract m_NStripsPhiCell * m_StripY / 2.0
  const double strip_localpos_y = m_StripY * strip_index + m_StripY / 2. - m_NStripsPhiCell * m_StripY / 2.0;

  CLHEP::Hep3Vector strip_localpos(m_StripXOffset, strip_localpos_y, strip_localpos_z);

  // Strip rotation
  const double phi = m_OffsetPhi + m_dPhi * segment_phi_bin;
  const double rotate = phi + m_OffsetRot;

  CLHEP::HepRotation rot;
  rot.rotateZ(rotate);
  strip_localpos = rot * strip_localpos;
  strip_localpos += ladder;

  location[0] = strip_localpos.x();
  location[1] = strip_localpos.y();
  location[2] = strip_localpos.z();

}

void CylinderGeomIntt::find_strip_index_values(const int segment_z_bin, const double yin, const double zin, int &strip_y_index, int &strip_z_index)
{
  // Given the location in y and z in sensor local coordinates, find the strip y and z index values

  // find the sensor type (inner or outer) from the segment_z_bin (location of sensor on ladder)
  const int itype = segment_z_bin % 2;
  if (itype != 0 && itype != 1)
  {
    cout << "Problem: itype = " << itype << endl;
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
  strip_z_index = (int) (zup / strip_z);

  // get the strip y index
  double yup = (double) nstrips_y_sensor * m_StripY / 2.0 + ypos;
  strip_y_index = (int) (yup / m_StripY);

  /*
  cout << "segment_z_bin " << segment_z_bin << " ypos " << ypos << " zpos " << zpos << " zup " << zup << " yup " << yup << endl;
  cout << "      -- itype " << itype << " strip_y " << m_StripY << " strip_z " << strip_z << " nstrips_z_sensor " << nstrips_z_sensor 
       << " nstrips_y_sensor " << nstrips_y_sensor << endl;
  cout << "      --  strip_z_index " << strip_z_index << " strip_y_index " << strip_y_index << endl;
  */
}

void CylinderGeomIntt::find_strip_center_localcoords(const int segment_z_bin, const int strip_y_index, const int strip_z_index, double location[])
{
  // find the sensor type (inner or outer) from the segment_z_bin (location of sensor on ladder)
  const int itype = segment_z_bin % 2;
  if (itype != 0 && itype != 1)
  {
    cout << "Problem: itype = " << itype << endl;
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
