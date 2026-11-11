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

// Derived compact map. Schema7 additionally retains usable corners of eligible
// noncertified leaves; classification remains unchanged. Reads schema5 as well.
class PHGarfieldCompactMap
{
 public:
  using Geometry = PHGarfieldAdaptiveMap::Geometry;
  using Bounds = PHGarfieldAdaptiveMap::Bounds;
  using LeafClass = PHGarfieldAdaptiveMap::LeafClass;
  using FallbackReason = PHGarfieldAdaptiveMap::FallbackReason;
  static constexpr int LegacySchemaVersion = 5;
  static constexpr int LegacyCollectionSchemaVersion = 6;
  static constexpr int SchemaVersion = 7;
  // Collection is an outer-file provenance schema. Runtime arrays use schema7
  // (schema5 for legacy collections); source_file_sha256 binds the
  // embedded TSV of source files rather than one merged training ROOT file.
  static constexpr int CollectionSchemaVersion = 8;
  static constexpr const char* CollectionSchemaName = "PHG4TpcGarfieldCompactTransportLUTCollection";
  static constexpr const char* CollectionSourceName = "PHG4TpcGarfieldFullTransportLUTShardCollection";
  static constexpr const char* CollectionSourceKind = "exact_shard_conversion_and_offset_merge_v1";
  static constexpr std::uint32_t Invalid = std::numeric_limits<std::uint32_t>::max();
  static constexpr const char* SchemaName = "PHG4TpcGarfieldCompactTransportLUT";
  static constexpr const char* LegacyLayout = "primitive_u32_nodes12_corners32_values24_source_cell_order_v1";
  static constexpr const char* Layout = "primitive_u32_nodes12_corners32_values24_usable_uncertified_v2";
  static constexpr const char* LegacyConversionContract = "exhaustive_source_bounds_routing_corner_bits_v1";
  static constexpr const char* ConversionContract = "exhaustive_source_bounds_routing_usable_corner_bits_v2";
  static const char* LayoutForSchema(int version);
  // Engineering ceilings only; callers must explicitly opt in above defaults.
  // All indexed counts remain below the reserved u32 Invalid sentinel.
  static constexpr std::uint64_t HardMaximumSourcePoints = 1500000000ULL;
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
    std::uint64_t max_source_points{600000000};
  };
  struct Node
  {
    // Internal: tag0..2 is split axis; first/second are children.
    // Leaf: tag low2bits=3, class bits2..4, reason bits5..7; upperbits0.
    // first is corner-record ID if available; otherwise Invalid; second=0.
    // Schema5 only retains Certified. Schema7 may retain safe class0/4 corners.
    std::uint32_t first{Invalid}, second{0}, tag{0};
  };
  using Corners = std::array<std::uint32_t, 8>;
  using Value = std::array<double, 3>;
  struct SourceBinding
  {
    std::string file_sha256;
    std::string configuration_fingerprint_fnv1a64;
    std::string content_fingerprint_fnv1a64;
    std::string canonical_policy;
    std::string converter_sha256;
  };
  struct ExpectedSource
  {
    std::string file_sha256;
    std::string configuration_fingerprint_fnv1a64;
    std::string content_fingerprint_fnv1a64;
  };
  struct ConversionProof
  {
    std::uint64_t source_points{0}, source_cells{0}, roots_checked{0};
    std::uint64_t bounds_checked{0}, routes_checked{0}, certified_leaves_checked{0}, corner_values_checked{0};
    std::uint64_t uncertified_leaves_checked{0};
    std::string sha256;
  };
  struct Data
  {
    int schema_version{SchemaVersion};
    Geometry geometry;
    SourceBinding source;
    ConversionProof proof;
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
  // converter, reader and callers that configure explicit resource ceilings.
  static bool ValidateLimits(const Limits& limits, std::ostream& error);
  static bool Validate(const Data& data, const Limits& limits, std::ostream& error, bool complete = true);
  // Initialization-time fallback requirement; packed tags stay private here.
  static std::uint64_t CountNonCertifiedLeaves(const Data& data);
  // Validate schema4 again, then exhaustively prove the compact conversion.
  static bool Convert(PHGarfieldAdaptiveMap::Data& source, const SourceBinding& binding, Data& result,
                      const Limits& limits, std::ostream& error, bool complete = true,
                      bool retain_uncertified_corners = true);
  static bool VerifyConversion(const PHGarfieldAdaptiveMap::Data& source, const Data& data,
                               std::ostream& error);
  static bool FindLeaf(const Data& data, unsigned int side, double phi, double r, double z, Lookup& result);
  static bool Interpolate(const Data& data, const Lookup& cell, double phi, double r, double z,
                          PHGarfieldForwardMap::NodeValue& result);
  static bool HasUsableCorners(const Data& data, const Lookup& cell);
  // Explicit relaxation never upgrades a leaf's class or promises guard_um.
  // certified is true only when the original strict Interpolate path succeeds.
  static bool InterpolateRelaxed(const Data& data, const Lookup& cell, double phi, double r, double z,
                                 PHGarfieldForwardMap::NodeValue& result, bool& certified);
  static std::string ContentFingerprint(const Data& data);
  static std::string ProofCanonical(const Data& data);
  static std::string Sha256(const std::string& text);
  static std::string FileSha256(const std::string& path);
  // ROOT directory is GarfieldFullTransport/Runtime. Caller supplies geometry
  // parsed from SourceMetadata and trusted hashes from its frozen manifest.
  // Caller MUST authenticate the complete runtime file SHA256 against the
  // external conversion manifest before opening/reading it. Internal proof and
  // ExpectedSource do not authenticate an artifact that was jointly rewritten.
  static bool Read(TDirectory* directory, const Geometry& geometry, const ExpectedSource& expected, Data& result,
                   const Limits& limits, std::ostream& error, bool complete = true);
  static bool Write(TDirectory* directory, const Data& data, const Limits& limits, std::ostream& error);
};
static_assert(sizeof(PHGarfieldCompactMap::Node) == 12);
static_assert(sizeof(PHGarfieldCompactMap::Corners) == 32);
static_assert(sizeof(PHGarfieldCompactMap::Value) == 24);
#endif
