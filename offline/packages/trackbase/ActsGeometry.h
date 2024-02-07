#ifndef TRACKBASE_ACTSGEOMETRY_H
#define TRACKBASE_ACTSGEOMETRY_H

#include "ActsSurfaceMaps.h"

#include <Acts/Definitions/Units.hpp>

class TrkrCluster;

class ActsGeometry
{
 public:
  ActsGeometry() = default;
  ~ActsGeometry() {}

  void setGeometry(ActsTrackingGeometry& tGeometry)
  {
    m_tGeometry = tGeometry;
  }

  void setSurfMaps(ActsSurfaceMaps& surfMaps)
  {
    m_surfMaps = surfMaps;
  }

  ActsTrackingGeometry& geometry()
  {
    return m_tGeometry;
  }
  ActsSurfaceMaps& maps()
  {
    return m_surfMaps;
  }

  void set_drift_velocity(double vd) { _drift_velocity = vd; }
  double get_drift_velocity() { return _drift_velocity; }
  void set_extended_readout_time(float time) { _extened_readout_time = time; }
  float get_extended_readout_time() { return _extened_readout_time; }

  Eigen::Matrix<float, 3, 1> getGlobalPositionF(
      TrkrDefs::cluskey key,
      TrkrCluster* cluster);

  Acts::Vector3 getGlobalPosition(
      TrkrDefs::cluskey key,
      TrkrCluster* cluster);

  Acts::Vector3 getGlobalPositionTpc(
      TrkrDefs::cluskey key,
      TrkrCluster* cluster);

  Surface get_tpc_surface_from_coords(
      TrkrDefs::hitsetkey hitsetkey,
      Acts::Vector3 world,
      TrkrDefs::subsurfkey& subsurfkey);

  Acts::Transform3 makeAffineTransform(Acts::Vector3 rotation, Acts::Vector3 translation);

  Acts::Vector2 getLocalCoords(TrkrDefs::cluskey key, TrkrCluster* cluster);

 private:
  ActsTrackingGeometry m_tGeometry;
  ActsSurfaceMaps m_surfMaps;

  float _extened_readout_time = 0;  // ns
  double _drift_velocity = 8.0e-3;  // cm/ns
};

#endif
