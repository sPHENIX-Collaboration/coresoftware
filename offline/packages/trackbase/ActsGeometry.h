#ifndef TRACKBASE_ACTSGEOMETRY_H
#define TRACKBASE_ACTSGEOMETRY_H

#include "ActsSurfaceMaps.h"

#include <Acts/Definitions/Units.hpp>

class TrkrCluster;

class ActsGeometry
{
 public:
  ActsGeometry() = default;
  ~ActsGeometry() = default;

  void setGeometry(const ActsTrackingGeometry& tGeometry)
  {
    m_tGeometry = tGeometry;
  }

  void setSurfMaps(const ActsSurfaceMaps& surfMaps)
  {
    m_surfMaps = surfMaps;
  }

  //! const accessor
  const ActsTrackingGeometry& geometry() const
  {
    return m_tGeometry;
  }

  //! mutable accessor
  ActsTrackingGeometry& geometry()
  {
    return m_tGeometry;
  }

  //! const accessor
  const ActsSurfaceMaps& maps() const
  {
    return m_surfMaps;
  }

  //! mutable accessor
  ActsSurfaceMaps& maps()
  {
    return m_surfMaps;
  }

  void set_drift_velocity(double vd) { _drift_velocity = vd; }

  void set_tpc_tzero(double tz) { _tpc_tzero = tz; }
  double get_tpc_tzero() const { return _tpc_tzero; }

  double get_drift_velocity() const { return _drift_velocity; }

  Acts::Vector3 getGlobalPosition(
      TrkrDefs::cluskey key,
      TrkrCluster* cluster) const;

  Acts::Vector3 getGlobalPositionTpc(
      TrkrDefs::cluskey key,
      TrkrCluster* cluster) const;

  Acts::Vector3 getGlobalPositionTpc(
      const TrkrDefs::hitsetkey& hitsetkey, const TrkrDefs::hitkey& hitkey, const float& phi, const float& rad,
      const float& clockPeriod) const;

  Surface get_tpc_surface_from_coords(
      TrkrDefs::hitsetkey hitsetkey,
      Acts::Vector3 world,
      TrkrDefs::subsurfkey& subsurfkey) const ;

  Acts::Transform3 makeAffineTransform(Acts::Vector3 rotation, Acts::Vector3 translation) const;

  Acts::Vector2 getLocalCoords(TrkrDefs::cluskey key, TrkrCluster* cluster) const;

 private:
  ActsTrackingGeometry m_tGeometry;
  ActsSurfaceMaps m_surfMaps;
  double _drift_velocity = 8.0e-3;  // cm/ns
  double _tpc_tzero = 0.0;  // ns
};

#endif
