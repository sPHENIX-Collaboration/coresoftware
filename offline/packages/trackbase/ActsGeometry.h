#ifndef TRACKBASE_ACTSGEOMETRY_H
#define TRACKBASE_ACTSGEOMETRY_H

#include "ActsSurfaceMaps.h"

#include <Acts/Definitions/Units.hpp>

class TrkrCluster;

class ActsGeometry {

 public:
  ActsGeometry() = default;
  ~ActsGeometry() {} 

  void setGeometry(ActsTrackingGeometry& tGeometry) 
    { m_tGeometry = tGeometry; }

  void setSurfMaps(ActsSurfaceMaps& surfMaps)
    { m_surfMaps = surfMaps; }
  
  ActsTrackingGeometry& geometry()
    { return m_tGeometry; }
  ActsSurfaceMaps& maps() 
    { return m_surfMaps; }

  void set_drift_velocity(double vd) {_drift_velocity = vd;}
  double get_drift_velocity() {return _drift_velocity;}

  Eigen::Matrix<float,3,1> getGlobalPositionF(
      TrkrDefs:: cluskey key,       
      TrkrCluster* cluster);

  Acts::Vector3 getGlobalPosition(
      TrkrDefs:: cluskey key,       
      TrkrCluster* cluster);

  Acts::Vector3 getGlobalPositionTpc(
      TrkrDefs:: cluskey key,       
      TrkrCluster* cluster);

  Surface get_tpc_surface_from_coords(
      TrkrDefs::hitsetkey hitsetkey,
      Acts::Vector3 world,
      TrkrDefs::subsurfkey& subsurfkey);

 private:
  ActsTrackingGeometry m_tGeometry;
  ActsSurfaceMaps m_surfMaps;

  double _drift_velocity = 8.0e-3;  // cm/ns
};

#endif
