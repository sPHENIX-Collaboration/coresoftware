#ifndef TRACKBASE_ACTSSOURCELINK_H
#define TRACKBASE_ACTSSOURCELINK_H

#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <ActsExamples/EventData/GeometryContainers.hpp>
#include <ActsExamples/EventData/Index.hpp>

#include "TrkrDefs.h"

#include <cassert>
#include <iostream>

/// A source link that stores just an index.
///
/// This is intentionally kept as barebones as possible. The source link
/// is just a reference and will be copied, moved around, etc. often.
/// Keeping it small and separate from the actual, potentially large,
/// measurement data should result in better overall performance.
/// Using an index instead of e.g. a pointer, means source link and
/// measurement are decoupled and the measurement represenation can be
/// easily changed without having to also change the source link.
class ActsSourceLink final : public Acts::SourceLink {
 public:
  /// Construct from geometry identifier and index.
  constexpr ActsSourceLink(Acts::GeometryIdentifier gid, ActsExamples::Index idx)
    : SourceLink(gid), m_index(idx), m_cluskey(0) {}
  constexpr ActsSourceLink(Acts::GeometryIdentifier gid, ActsExamples::Index idx, TrkrDefs::cluskey cluskey)
    : SourceLink(gid), m_index(idx), m_cluskey(cluskey) {}

  // Construct an invalid source link. Must be default constructible to
  /// satisfy SourceLinkConcept.
  ActsSourceLink() : SourceLink{Acts::GeometryIdentifier{}} {m_cluskey = 0;}
  ActsSourceLink(const ActsSourceLink&) = default;
  ActsSourceLink(ActsSourceLink&&) = default;
  ActsSourceLink& operator=(const ActsSourceLink&) = default;
  ActsSourceLink& operator=(ActsSourceLink&&) = default;

  /// Access the index.
  constexpr ActsExamples::Index index() const { return m_index; }
  constexpr TrkrDefs::cluskey cluskey() const { return m_cluskey; }
 private:
  ActsExamples::Index m_index;
  TrkrDefs::cluskey m_cluskey;

  friend constexpr bool operator==(const ActsSourceLink& lhs,
                                   const ActsSourceLink& rhs) {
    return (lhs.geometryId() == rhs.geometryId()) and
           (lhs.m_index == rhs.m_index) and
           (lhs.m_cluskey == rhs.m_cluskey);
  }
  friend constexpr bool operator!=(const ActsSourceLink& lhs,
                                   const ActsSourceLink& rhs) {
    return not(lhs == rhs);
  }
};

/// Container of index source links.
///
/// Since the source links provide a `.geometryId()` accessor, they can be
/// stored in an ordered geometry container.
using ActsSourceLinkContainer =
  ActsExamples::GeometryIdMultiset<std::reference_wrapper<const ActsSourceLink>>;
/// Accessor for the above source link container
///
/// It wraps up a few lookup methods to be used in the Combinatorial Kalman
/// Filter
using ActsSourceLinkAccessor =
  ActsExamples::GeometryIdMultisetAccessor<std::reference_wrapper<const ActsSourceLink>>;

#endif
