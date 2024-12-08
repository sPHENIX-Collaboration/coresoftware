#include "CylinderGeomInttHelper.h"

#include <trackbase/ActsGeometry.h>          // for ActsGeometry
#include <trackbase/ActsTrackingGeometry.h>  // for ActsTrackingGeometry

#include <Acts/Definitions/Units.hpp>

#include <CLHEP/Vector/Rotation.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <algorithm>
#include <cmath>
#include <memory>  // for __shared_ptr_access
#include <utility>

TVector3 CylinderGeomInttHelper::get_world_from_local_coords(const Surface& surface, ActsGeometry* tGeometry, const TVector3& local)
{
  Acts::Vector3 loc(local.x(), local.y(), local.z());
  loc *= Acts::UnitConstants::cm;

  Acts::Vector3 glob = surface->transform(tGeometry->geometry().getGeoContext()) * loc;
  glob /= Acts::UnitConstants::cm;
  return TVector3(glob(0), glob(1), glob(2));
}

TVector3 CylinderGeomInttHelper::get_world_from_local_coords(const Surface& surface, ActsGeometry* tGeometry, const TVector2& local)
{
  Acts::Vector2 actslocal;
  actslocal(0) = local.X();
  actslocal(1) = local.Y();
  actslocal *= Acts::UnitConstants::cm;

  /// Acts requires a dummy vector to be passed in the arg list
  auto global = surface->localToGlobal(tGeometry->geometry().getGeoContext(),
                                       actslocal,
                                       Acts::Vector3(1, 1, 1));

  global /= Acts::UnitConstants::cm;

  TVector3 ret;
  ret[0] = global(0);
  ret[1] = global(1);
  ret[2] = global(2);

  return ret;
}

TVector3 CylinderGeomInttHelper::get_local_from_world_coords(const Surface& surface, ActsGeometry* tGeometry, TVector3 world)
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

void CylinderGeomInttHelper::find_segment_center(const Surface& surface, ActsGeometry* tGeometry, double location[])
{
  TVector2 local(0.0, 0.0);

  TVector3 global = get_world_from_local_coords(surface, tGeometry, local);
  location[0] = global.X();
  location[1] = global.Y();
  location[2] = global.Z();
  return;
}

void
CylinderGeomInttHelper::find_strip_center (
  const Surface& surface,
  ActsGeometry* tGeometry,
  const int segment_z_bin,
  const int segment_phi_bin,
  const int strip_column,
  const int strip_index,
  double* location,
  CylinderGeomIntt& fren
) {
  // Ladder
  find_segment_center(surface, tGeometry, location);
  CLHEP::Hep3Vector ladder(location[0], location[1], location[2]);

  // Strip
  const int itype = segment_z_bin % 2;
  const double strip_z = fren.m_StripZ[itype];
  const int nstrips_z_sensor = fren.m_NStripsZSensor[itype];

  const double strip_localpos_z = strip_z * (strip_column % nstrips_z_sensor) - strip_z / 2. * nstrips_z_sensor + strip_z / 2.;
  // distance from bottom of sensor = m_StripY*strip_index +m_StripY/2.0, then subtract m_NStripsPhiCell * m_StripY / 2.0
  const double strip_localpos_y = fren.m_StripY * strip_index + fren.m_StripY / 2. - fren.m_NStripsPhiCell * fren.m_StripY / 2.0;

  CLHEP::Hep3Vector strip_localpos(fren.m_StripXOffset, strip_localpos_y, strip_localpos_z);

  // Strip rotation
  const double phi = fren.m_OffsetPhi + fren.m_dPhi * segment_phi_bin;
  const double rotate = phi + fren.m_OffsetRot;

  CLHEP::HepRotation rot;
  rot.rotateZ(rotate);
  strip_localpos = rot * strip_localpos;
  strip_localpos += ladder;

  location[0] = strip_localpos.x();
  location[1] = strip_localpos.y();
  location[2] = strip_localpos.z();
}

