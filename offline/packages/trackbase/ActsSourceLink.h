#ifndef TRACKBASE_ACTSSOURCELINK_H
#define TRACKBASE_ACTSSOURCELINK_H

#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include "TrkrDefs.h"

#include <cassert>
#include <iostream>

/// A source link that stores just an index.
/// Using an index instead of e.g. a pointer, means source link and
/// measurement are decoupled and the measurement represenation can be
/// easily changed without having to also change the source link.
class ActsSourceLink final
{
 public:
  using Index = uint8_t;

  /// Construct from geometry identifier and index.
  constexpr ActsSourceLink(Acts::GeometryIdentifier gid, Index idx)
    : m_geometryId(gid)
    , m_index(idx)
    , m_cluskey(0)
  {
  }
  constexpr ActsSourceLink(Acts::GeometryIdentifier gid, Index idx, TrkrDefs::cluskey cluskey)
    : m_geometryId(gid)
    , m_index(idx)
    , m_cluskey(cluskey)
  {
  }

  // Construct an invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  ActsSourceLink()
    : m_geometryId{Acts::GeometryIdentifier{}}
    , m_index(UINT8_MAX)
    , m_cluskey(UINT64_MAX)
  {
  }

  ActsSourceLink(const ActsSourceLink&) = default;
  ActsSourceLink(ActsSourceLink&&) = default;
  ActsSourceLink& operator=(const ActsSourceLink&) = default;
  ActsSourceLink& operator=(ActsSourceLink&&) = default;

  /// Access the index.
  constexpr Index index() const { return m_index; }
  constexpr TrkrDefs::cluskey cluskey() const { return m_cluskey; }
  constexpr Acts::GeometryIdentifier geometryId() const { return m_geometryId; }

  struct SurfaceAccessor
  {
    const Acts::TrackingGeometry& trackingGeometry;
    const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const
    {
      const auto& sl = sourceLink.get<ActsSourceLink>();
      return trackingGeometry.findSurface(sl.geometryId());
    }
  };

 private:
  Acts::GeometryIdentifier m_geometryId;
  Index m_index;
  TrkrDefs::cluskey m_cluskey;

  friend constexpr bool operator==(const ActsSourceLink& lhs,
                                   const ActsSourceLink& rhs)
  {
    return (lhs.geometryId() == rhs.geometryId()) and
           (lhs.m_index == rhs.m_index) and
           (lhs.m_cluskey == rhs.m_cluskey);
  }
  friend constexpr bool operator!=(const ActsSourceLink& lhs,
                                   const ActsSourceLink& rhs)
  {
    return not(lhs == rhs);
  }
};

#endif
