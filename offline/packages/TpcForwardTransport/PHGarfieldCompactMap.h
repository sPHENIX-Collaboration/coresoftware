#ifndef PHGARFIELD_COMPACTMAP_H
#define PHGARFIELD_COMPACTMAP_H

#include <array>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>
#include "PHGarfieldAdaptiveMap.h"

class TDirectory;

// Adaptive transport LUT: runtime routing, retained corner values and policy.
class PHGarfieldCompactMap
{
 public:
  using Geometry = PHGarfieldAdaptiveMap::Geometry;
  using Bounds = PHGarfieldAdaptiveMap::Bounds;
  using LeafClass = PHGarfieldAdaptiveMap::LeafClass;
  using FallbackReason = PHGarfieldAdaptiveMap::FallbackReason;
  static constexpr int SchemaVersion = 9;
  static constexpr std::uint32_t Invalid = std::numeric_limits<std::uint32_t>::max();
  static constexpr const char* SchemaName = "PHG4TpcGarfieldAdaptiveTransportLUT";
  static constexpr const char* Layout = "primitive_u32_nodes12_corners32_values24_usable_uncertified_v2";
  // Engineering ceilings only; callers must explicitly opt in above defaults.
  // All indexed counts remain below the reserved u32 Invalid sentinel.
  static constexpr std::uint64_t HardMaximumNodes = 250000000ULL;
  static constexpr std::uint64_t HardMaximumValues = 1500000000ULL;
  static constexpr std::uint64_t HardMaximumCorners = 150000000ULL;
  static constexpr std::uint64_t HardMaximumRoots = 1000000ULL;
  static constexpr std::uint64_t HardMaximumBytes = 128ULL * 1024 * 1024 * 1024;

  struct Limits
  {
    std::uint64_t max_nodes{100000000};
    std::uint64_t max_values{600000000};
    std::uint64_t max_corners{50000000};
    std::uint64_t max_roots{1000000};
    std::uint64_t max_bytes{32ULL * 1024 * 1024 * 1024};
  };
  struct Node
  {
    // Internal: tag0..2 is split axis; first/second are children.
    // Leaf: tag low2bits=3, class bits2..4, reason bits5..7; upperbits0.
    // first is corner-record ID if available; otherwise Invalid; second=0.
    // Certified and eligible noncertified leaves retain safe corner values.
    std::uint32_t first{Invalid}, second{0}, tag{0};
  };
  using Corners = std::array<std::uint32_t, 8>;
  using Value = std::array<double, 3>;
  struct Data
  {
    Geometry geometry;
    std::string configuration_fingerprint_fnv1a64;
    std::string canonical_policy;
    unsigned int max_depth{0}, max_nodes_per_root{0}, max_queries_per_root{0};
    double guard_um{0.}, maximum_drift_time_ns{0.};
    std::uint64_t first_root{0}, end_root{0};
    std::vector<Node> nodes;
    std::vector<Corners> corners;
    std::vector<Value> values;
    std::vector<std::uint32_t> root_lookup;
  };
  struct Lookup
  {
    std::uint32_t cell_id{Invalid}, corner_record{Invalid};
    unsigned int depth{0};
    int leaf_class{-1}, fallback_reason{-1};
    Bounds bounds;
  };

  static bool EstimatedBytes(std::uint64_t nodes, std::uint64_t corners, std::uint64_t values,
                             std::uint64_t roots, std::uint64_t& bytes);
  // Scalar checks only: no allocation, ROOT access or traversal. Shared by the
  // packing tools, reader and callers configuring explicit resource ceilings.
  static bool ValidateLimits(const Limits& limits, std::ostream& error);
  static bool Validate(const Data& data, const Limits& limits, std::ostream& error, bool complete = true);
  // Initialization-time fallback requirement; packed tags stay private here.
  static std::uint64_t CountNonCertifiedLeaves(const Data& data);
  // Construction tools only: pack the in-memory builder output, or append a
  // compatible contiguous shard. No intermediate ROOT format is required.
  static bool Pack(const PHGarfieldAdaptiveMap::Data& source, const std::string& configuration_fingerprint,
                   Data& result, const Limits& limits, std::ostream& error, bool complete = true);
  static bool Append(Data& result, const Data& shard, const Limits& limits, std::ostream& error);
  static bool FindLeaf(const Data& data, unsigned int side, double phi, double r, double z, Lookup& result);
  static bool Interpolate(const Data& data, const Lookup& cell, double phi, double r, double z,
                          PHGarfieldForwardMap::NodeValue& result);
  static bool HasUsableCorners(const Data& data, const Lookup& cell);
  // Explicit relaxation never upgrades a leaf's class or promises guard_um.
  // certified is true only when the original strict Interpolate path succeeds.
  static bool InterpolateRelaxed(const Data& data, const Lookup& cell, double phi, double r, double z,
                                 PHGarfieldForwardMap::NodeValue& result, bool& certified);
  static std::string ContentFingerprint(const Data& data);
  // ROOT directory is GarfieldFullTransport/Runtime; geometry is parsed from
  // its parent Metadata. The content fingerprint binds geometry and arrays.
  static bool Read(TDirectory* directory, const Geometry& geometry, Data& result,
                   const Limits& limits, std::ostream& error, bool complete = true);
  static bool Write(TDirectory* directory, const Data& data, const Limits& limits, std::ostream& error);
};
static_assert(sizeof(PHGarfieldCompactMap::Node) == 12);
static_assert(sizeof(PHGarfieldCompactMap::Corners) == 32);
static_assert(sizeof(PHGarfieldCompactMap::Value) == 24);
#endif
