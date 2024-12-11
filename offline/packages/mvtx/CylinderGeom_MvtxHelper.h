#ifndef MVTX_CYLINDERGEOMMVTXHELPER_H
#define MVTX_CYLINDERGEOMMVTXHELPER_H

#include <g4detectors/PHG4CylinderGeom.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <TVector3.h>

#include <iostream>

class CylinderGeom_MvtxHelper
{
 public:

  // our own - no override
  TVector3 static
  get_local_from_world_coords (
    Surface const& surface,
    ActsGeometry* tGeometry,
    TVector3 world
  );

  TVector3 static
  get_world_from_local_coords (
    Surface const& surface,
    ActsGeometry* tGeometry,
    TVector2 const& local
  );

  TVector3 static
  get_world_from_local_coords (
    Surface const& surface,
    ActsGeometry* tGeometry,
    TVector3 const& local
  );

  void static
  find_sensor_center (
    Surface const& surface,
    ActsGeometry* tGeometry,
    double* location
  );

};

#endif
