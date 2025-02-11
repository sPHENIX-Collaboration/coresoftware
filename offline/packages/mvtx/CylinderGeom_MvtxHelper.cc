#include "CylinderGeom_MvtxHelper.h"

#include "SegmentationAlpide.h"

#include <Acts/Definitions/Units.hpp>

#include <TRotation.h>
#include <TVector3.h>

#include <cmath>
#include <ostream>  // for operator<<, basic_ostream::operator<<, basic_...

TVector3
CylinderGeom_MvtxHelper::get_local_from_world_coords (
  const Surface& surface,
  ActsGeometry* tGeometry,
  TVector3 world
) {
  Acts::Vector3 global;
  global(0) = world[0];
  global(1) = world[1];
  global(2) = world[2];

  global *= Acts::UnitConstants::cm;

  Acts::Vector3 local = surface->transform(tGeometry->geometry().getGeoContext()).inverse() * global;
  local /= Acts::UnitConstants::cm;

  /// The Acts transform swaps a few of the coordinates
  return TVector3(local(0), local(2) * -1, local(1));
}

TVector3
CylinderGeom_MvtxHelper::get_world_from_local_coords (
  const Surface& surface,
  ActsGeometry* tGeometry,
  const TVector2& local
) {
  Acts::Vector2 actslocal;
  actslocal(0) = local.X();
  actslocal(1) = local.Y();
  actslocal *= Acts::UnitConstants::cm;

  Acts::Vector3 global;
  /// Acts requires a dummy vector to be passed in the arg list
  global = surface->localToGlobal (
    tGeometry->geometry().getGeoContext(),
    actslocal,
    Acts::Vector3(1, 1, 1)
  );
  global /= Acts::UnitConstants::cm;

  TVector3 res;
  res[0] = global(0);
  res[1] = global(1);
  res[2] = global(2);

  return res;
}

TVector3
CylinderGeom_MvtxHelper::get_world_from_local_coords (
  const Surface& surface,
  ActsGeometry* tGeometry,
  const TVector3& local
) {
  Acts::Vector3 loc(local.x(), local.y(), local.z());
  loc *= Acts::UnitConstants::cm;

  Acts::Vector3 glob = surface->transform(tGeometry->geometry().getGeoContext()) * loc;
  glob /= Acts::UnitConstants::cm;

  return TVector3(glob(0), glob(1), glob(2));
}

void CylinderGeom_MvtxHelper::find_sensor_center (
		const Surface& surface,
		ActsGeometry* tGeometry,
		double* location
) {
  TVector2 sensor_local(0.0, 0.0);

  TVector3 sensor_world = get_world_from_local_coords(surface, tGeometry, sensor_local);

  location[0] = sensor_world.X();
  location[1] = sensor_world.Y();
  location[2] = sensor_world.Z();

  return;
}

