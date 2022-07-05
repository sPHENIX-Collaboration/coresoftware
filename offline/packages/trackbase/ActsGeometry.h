#ifndef TRACKBASE_ACTSGEOMETRY_H
#define TRACKBASE_ACTSGEOMETRY_H

#include "ActsSurfaceMaps.h"
#include <Acts/Definitions/Units.hpp>

class TrkrCluster;

class ActsGeometry {

 public:
  ActsGeometry() = default;
  ~ActsGeometry() {} 

  void setGeometry(const ActsTrackingGeometry tGeometry) 
    { m_tGeometry = tGeometry; }

  void setSurfMaps(const ActsSurfaceMaps surfMaps)
    { m_surfMaps = surfMaps; }
  
  ActsTrackingGeometry geometry() const 
    { return m_tGeometry; }
  ActsSurfaceMaps maps() const 
    { return m_surfMaps; }

  Eigen::Matrix<float,3,1> getGlobalPositionF(
      TrkrDefs:: cluskey key,       
      TrkrCluster* cluster) const;

  Acts::Vector3 getGlobalPosition(
      TrkrDefs:: cluskey key,       
      TrkrCluster* cluster) const;

  Surface get_tpc_surface_from_coords(
      TrkrDefs::hitsetkey hitsetkey,
      Acts::Vector3 world,
      TrkrDefs::subsurfkey& subsurfkey) const;

 private:
  ActsTrackingGeometry m_tGeometry;
  ActsSurfaceMaps m_surfMaps;

};

#endif
