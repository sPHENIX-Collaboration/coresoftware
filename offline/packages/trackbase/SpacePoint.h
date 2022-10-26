#ifndef TRACKBASE_SPACEPOINT_H
#define TRACKBASE_SPACEPOINT_H

#include <memory>
#include "trackbase/TrkrDefs.h"

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Seeding/Seed.hpp>

/**
 * A struct for Acts to take cluster information for seeding
 */
struct SpacePoint {
  TrkrDefs::cluskey m_clusKey;
  double m_x;
  double m_y;
  double m_z;
  double m_r;
  Acts::GeometryIdentifier m_geoId;
  double m_varianceR;
  double m_varianceZ;
  
  TrkrDefs::cluskey Id() const { return m_clusKey; }

  /// These are needed by Acts
  double x() const { return m_x; }
  double y() const { return m_y; }
  double z() const { return m_z; }
  double r() const { return m_r; }
  double varianceR() const { return m_varianceR; }
  double varianceZ() const { return m_varianceZ; }

};

/// This is needed by the Acts seedfinder 
inline bool operator==(SpacePoint a, SpacePoint b) {
  return (a.m_clusKey == b.m_clusKey);
}

using SpacePointPtr = std::unique_ptr<SpacePoint>;
using SeedContainer = std::vector<Acts::Seed<SpacePoint>>;

#endif
