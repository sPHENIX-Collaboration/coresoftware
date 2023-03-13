#ifndef TRACKBASE_ACTSSOURCELINK_H
#define TRACKBASE_ACTSSOURCELINK_H

#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include "TrkrDefs.h"

#include <cassert>
#include <iostream>

/// A source link that stores just an index.
/// Using an index instead of e.g. a pointer, means source link and
/// measurement are decoupled and the measurement represenation can be
/// easily changed without having to also change the source link.
class ActsSourceLink final : public Acts::SourceLink
{
 public:
  using Index = uint8_t;

  /// Construct from geometry identifier and index.
  constexpr ActsSourceLink(Acts::GeometryIdentifier gid, Index idx)
    : SourceLink(gid)
    , m_index(idx)
    , m_cluskey(0)
  {
  }
  constexpr ActsSourceLink(Acts::GeometryIdentifier gid, Index idx, TrkrDefs::cluskey cluskey)
    : SourceLink(gid)
    , m_index(idx)
    , m_cluskey(cluskey)
  {
  }

  // Construct an invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  ActsSourceLink()
    : SourceLink{Acts::GeometryIdentifier{}}
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

 private:
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
