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

// Derived schema5. The complete schema4 training map remains immutable evidence.
// Only certified corner values survive; row index preserves the training cell ID.
class PHGarfieldCompactMap
{
 public:
  using Geometry = PHGarfieldAdaptiveMap::Geometry;
  using Bounds = PHGarfieldAdaptiveMap::Bounds;
  using LeafClass = PHGarfieldAdaptiveMap::LeafClass;
  using FallbackReason = PHGarfieldAdaptiveMap::FallbackReason;
  static constexpr int SchemaVersion = 5;
  // Collection is an outer-file provenance schema. Runtime arrays still use
  // schema5 and the identical reader/interpolator; source_file_sha256 binds the
  // embedded TSV of source files rather than one merged training ROOT file.
  static constexpr int CollectionSchemaVersion = 6;
  static constexpr const char* CollectionSchemaName = "PHG4TpcGarfieldCompactTransportLUTCollection";
  static constexpr const char* CollectionSourceName = "PHG4TpcGarfieldFullTransportLUTShardCollection";
  static constexpr const char* CollectionSourceKind = "exact_shard_conversion_and_offset_merge_v1";
  static constexpr std::uint32_t Invalid = std::numeric_limits<std::uint32_t>::max();
  static constexpr const char* SchemaName = "PHG4TpcGarfieldCompactTransportLUT";
  static constexpr const char* Layout = "primitive_u32_nodes12_corners32_values24_source_cell_order_v1";
  static constexpr const char* ConversionContract = "exhaustive_source_bounds_routing_corner_bits_v1";
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
    // Certified first is corner-record ID; direct first=Invalid; second=0.
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
    std::string sha256;
  };
  struct Data
  {
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
  static bool ValidateLimits(const Limits&, std::ostream&);
  static bool Validate(const Data&, const Limits&, std::ostream&, bool require_complete = true);
  // Initialization-time fallback requirement; packed tags stay private here.
  static std::uint64_t CountNonCertifiedLeaves(const Data&);
  // Validate schema4 again, then exhaustively prove the compact conversion.
  static bool Convert(PHGarfieldAdaptiveMap::Data&, const SourceBinding&, Data&,
                      const Limits&, std::ostream&, bool require_complete = true);
  static bool VerifyConversion(const PHGarfieldAdaptiveMap::Data&, const Data&,
                               std::ostream&);
  static bool FindLeaf(const Data&, unsigned int side, double phi, double r, double z, Lookup&);
  static bool Interpolate(const Data&, const Lookup&, double phi, double r, double z,
                          PHGarfieldForwardMap::NodeValue&);
  static std::string ContentFingerprint(const Data&);
  static std::string ProofCanonical(const Data&);
  static std::string Sha256(const std::string&);
  static std::string FileSha256(const std::string&);
  // ROOT directory is GarfieldFullTransport/Runtime. Caller supplies geometry
  // parsed from SourceMetadata and trusted hashes from its frozen manifest.
  // Caller MUST authenticate the complete runtime file SHA256 against the
  // external conversion manifest before opening/reading it. Internal proof and
  // ExpectedSource do not authenticate an artifact that was jointly rewritten.
  static bool Read(TDirectory*, const Geometry&, const ExpectedSource&, Data&,
                   const Limits&, std::ostream&, bool require_complete = true);
  static bool Write(TDirectory*, const Data&, const Limits&, std::ostream&);
};
static_assert(sizeof(PHGarfieldCompactMap::Node) == 12);
static_assert(sizeof(PHGarfieldCompactMap::Corners) == 32);
static_assert(sizeof(PHGarfieldCompactMap::Value) == 24);
#endif
