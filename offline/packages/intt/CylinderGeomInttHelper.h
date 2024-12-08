#ifndef INTT_CYLINDERGEOMINTTHELPER_H
#define INTT_CYLINDERGEOMINTTHELPER_H

#include "CylinderGeomIntt.h"

#include <g4detectors/PHG4CylinderGeom.h>

#include <trackbase/ActsGeometry.h>

#include <TVector2.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>

class CylinderGeomInttHelper
{
 public:
  CylinderGeomInttHelper() = default;

  // our own
  void static
  find_segment_center (
    const Surface& surface,
    ActsGeometry* tGeometry,
    double location[]
  );

  void static
  find_strip_center (
    Surface const& surface,
    ActsGeometry* tGeometry,
    const int segment_z_bin,
    const int segment_phi_bin,
    const int strip_column,
    const int strip_index,
    double* location,
	CylinderGeomIntt& fren
  );
 
  TVector3 static
  get_world_from_local_coords (
    const Surface& surface,
    ActsGeometry* tGeometry,
    const TVector2& local
  );

  TVector3 static
  get_world_from_local_coords (
    const Surface& surface,
    ActsGeometry* tGeometry,
    const TVector3& local
  );

  TVector3 static
  get_local_from_world_coords (
    const Surface& surface,
    ActsGeometry* tGeometry,
    TVector3 world
  );
};

#endif
