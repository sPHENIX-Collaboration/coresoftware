#include "PHGarfieldCompactMap.h"
#include <TBasket.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <openssl/evp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <locale>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include "PHGarfieldCompactMapInternal.h"
#include "PHGarfieldFingerprint.h"
#include "PHGarfieldPrimitiveTree.h"

namespace
{
  using Map = PHGarfieldCompactMap;
  using Training = PHGarfieldAdaptiveMap;
  constexpr double TwoPi = 2. * std::numbers::pi;
  using PHGarfieldCompactMapInternal::binding_valid;
  using PHGarfieldCompactMapInternal::class_reason;
  using PHGarfieldCompactMapInternal::hex;
  using PHGarfieldCompactMapInternal::kind;
  using PHGarfieldCompactMapInternal::leaf_class;
  using PHGarfieldCompactMapInternal::reason;
  using PHGarfieldCompactMapInternal::require;
  class Digest
  {
   public:
    Digest()
    {
      if (!ctx || EVP_DigestInit_ex(ctx.get(), EVP_sha256(), nullptr) != 1)
      {
        throw std::runtime_error("SHA256 initialization failed");
      }
    }
    void add(const char* p, std::size_t n) const
    {
      if (EVP_DigestUpdate(ctx.get(), p, n) != 1)
      {
        throw std::runtime_error("SHA256 update failed");
      }
    }
    std::string finish() const
    {
      std::array<unsigned char, EVP_MAX_MD_SIZE> bytes{};
      unsigned int n = 0;
      if (EVP_DigestFinal_ex(ctx.get(), bytes.data(), &n) != 1 || n != 32)
      {
        throw std::runtime_error("SHA256 finalize failed");
      }
      std::ostringstream out;
      out << std::hex << std::setfill('0');
      for (unsigned int i = 0; i < n; ++i)
      {
        out << std::setw(2) << static_cast<unsigned int>(bytes[i]);
      }
      return out.str();
    }

   private:
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> ctx{EVP_MD_CTX_new(), EVP_MD_CTX_free};
  };
}  // namespace

const char* PHGarfieldCompactMap::LayoutForSchema(int version)
{
  if (version == LegacySchemaVersion || version == LegacyCollectionSchemaVersion)
  {
    return LegacyLayout;
  }
  return version == SchemaVersion || version == CollectionSchemaVersion ? Layout : nullptr;
}

std::string PHGarfieldCompactMap::Sha256(const std::string& text)
{
  Digest digest;
  digest.add(text.data(), text.size());
  return digest.finish();
}
std::string PHGarfieldCompactMap::FileSha256(const std::string& path)
{
  std::ifstream input(path, std::ios::binary);
  if (!input)
  {
    throw std::runtime_error("cannot hash " + path);
  }
  Digest digest;
  std::array<char, 1024 * 1024> block{};
  while (input)
  {
    input.read(block.data(), block.size());
    const auto n = input.gcount();
    if (n > 0)
    {
      digest.add(block.data(), n);
    }
  }
  if (!input.eof())
  {
    throw std::runtime_error("SHA256 file read failed " + path);
  }
  return digest.finish();
}
bool PHGarfieldCompactMap::EstimatedBytes(std::uint64_t nodes, std::uint64_t corners, std::uint64_t values, std::uint64_t roots, std::uint64_t& bytes)
{
  bytes = 1024;
  // Conservative growth+validation allowance, above exact12/32/24/4-byte arrays.
  for (const auto& item : {std::make_pair(nodes, 32ULL), std::make_pair(corners, 96ULL), std::make_pair(values, 64ULL), std::make_pair(roots, 16ULL)})
  {
    if (item.first > (std::numeric_limits<std::uint64_t>::max() - bytes) / item.second)
    {
      return false;
    }
    bytes += item.first * item.second;
  }
  return true;
}

bool PHGarfieldCompactMap::ValidateLimits(const Limits& limits, std::ostream& error)
{
  return require(limits.max_source_points && limits.max_source_points <= HardMaximumSourcePoints &&
                     limits.max_nodes && limits.max_nodes <= HardMaximumNodes &&
                     limits.max_values && limits.max_values <= HardMaximumValues &&
                     limits.max_corners && limits.max_corners <= HardMaximumCorners &&
                     limits.max_roots && limits.max_roots <= HardMaximumRoots &&
                     limits.max_bytes >= 1024 && limits.max_bytes <= HardMaximumBytes,
                 error, "invalid or excessive compact resource limits");
}

std::uint64_t PHGarfieldCompactMap::CountNonCertifiedLeaves(const Data& data)
{
  std::uint64_t count = 0;
  for (const auto& node : data.nodes)
  {
    if (kind(node) == 3U && leaf_class(node) != static_cast<int>(LeafClass::Certified))
    {
      ++count;
    }
  }
  return count;
}

bool PHGarfieldCompactMap::Validate(const Data& data, const Limits& limits, std::ostream& error, bool complete)
{
  if (!ValidateLimits(limits, error))
  {
    return false;
  }
  const auto roots = Training::RootCount(data.geometry);
  std::uint64_t bytes = 0;
  if (!require((data.schema_version == LegacySchemaVersion || data.schema_version == SchemaVersion) &&
                   roots && roots <= limits.max_roots && binding_valid(data.source) &&
                   data.first_root <= data.end_root && data.end_root <= roots && (!complete || (!data.first_root && data.end_root == roots)) &&
                   data.root_lookup.size() == roots && data.nodes.size() <= limits.max_nodes && data.corners.size() <= limits.max_corners && data.values.size() <= limits.max_values &&
                   data.nodes.capacity() <= limits.max_nodes && data.corners.capacity() <= limits.max_corners && data.values.capacity() <= limits.max_values && data.root_lookup.capacity() <= limits.max_roots &&
                   EstimatedBytes(data.nodes.capacity(), data.corners.capacity(), data.values.capacity(), data.root_lookup.capacity(), bytes) && bytes <= limits.max_bytes &&
                   data.max_depth <= 63 && data.max_nodes_per_root && data.max_queries_per_root && std::isfinite(data.guard_um) && data.guard_um > 0. &&
                   std::isfinite(data.maximum_drift_time_ns) && data.maximum_drift_time_ns >= 0.,
               error, "invalid geometry/binding/count/capacity/budget"))
  {
    return false;
  }
  const auto& proof = data.proof;
  if (!require(proof.source_points >= data.values.size() && proof.source_points <= limits.max_source_points && proof.source_cells == data.nodes.size() &&
                   proof.roots_checked == data.end_root - data.first_root && proof.bounds_checked == data.nodes.size() && proof.routes_checked == data.nodes.size() &&
                   proof.certified_leaves_checked + proof.uncertified_leaves_checked == data.corners.size() &&
                   (data.schema_version != LegacySchemaVersion || proof.uncertified_leaves_checked == 0) &&
                   proof.corner_values_checked == 8ULL * data.corners.size() &&
                   hex(proof.sha256, 64) && proof.sha256 == Sha256(ProofCanonical(data)),
               error, "missing/inconsistent conversion proof"))
  {
    return false;
  }
  for (const auto& value : data.values)
  {
    if (!require(std::all_of(value.begin(), value.end(), [](double v)
                             { return std::isfinite(v); }) &&
                     value[2] >= 0. && value[2] <= data.maximum_drift_time_ns,
                 error, "invalid compact endpoint value"))
    {
      return false;
    }
  }
  std::vector<unsigned char> visited(data.nodes.size(), 0);
  std::vector<unsigned char> used_corners(data.corners.size(), 0);
  std::vector<unsigned char> used_values(data.values.size(), 0);
  std::uint64_t visited_count = 0;
  std::uint64_t certified_count = 0;
  std::uint64_t uncertified_count = 0;
  std::function<bool(std::uint32_t, unsigned int, const Bounds&)> visit;
  visit = [&](std::uint32_t id, unsigned int depth, const Bounds& bounds)
  {
    if (!require(id < data.nodes.size() && !visited[id] && depth <= data.max_depth, error, "missing/duplicate/cyclic compact node"))
    {
      return false;
    }
    visited[id] = 1;
    ++visited_count;
    const auto& n = data.nodes[id];
    for (unsigned int a = 0; a < 3; ++a)
    {
      if (!require(bounds.hi[a] > bounds.lo[a], error, "empty compact cell extent"))
      {
        return false;
      }
    }
    if (kind(n) < 3)
    {
      const unsigned int axis = kind(n);
      const double mid = .5 * (bounds.lo[axis] + bounds.hi[axis]);
      if (!require(n.tag == axis && depth < data.max_depth && mid > bounds.lo[axis] && mid < bounds.hi[axis], error, "invalid compact split/reserved bits"))
      {
        return false;
      }
      return visit(n.first, depth + 1, Training::ChildBounds(bounds, axis, 0)) && visit(n.second, depth + 1, Training::ChildBounds(bounds, axis, 1));
    }
    const int c = leaf_class(n);
    const int r = reason(n);
    if (!require(n.tag <= 255 && class_reason(c, r) && n.second == 0, error, "invalid compact leaf/reason/reserved bits"))
    {
      return false;
    }
    if (c != static_cast<int>(LeafClass::Certified) && n.first == Invalid)
    {
      return true;
    }
    if (c != static_cast<int>(LeafClass::Certified) &&
        !require(data.schema_version == SchemaVersion &&
                     (c == static_cast<int>(LeafClass::Uncertified) || c == static_cast<int>(LeafClass::NumericallyUnresolved)) &&
                     r != static_cast<int>(FallbackReason::AngularBranch),
                 error, "ineligible or legacy noncertified leaf has interpolation payload"))
    {
      return false;
    }
    if (!require(n.first < data.corners.size() && !used_corners[n.first], error, "invalid/duplicate corner record"))
    {
      return false;
    }
    used_corners[n.first] = 1;
    if (c == static_cast<int>(LeafClass::Certified))
    {
      ++certified_count;
    }
    else
    {
      ++uncertified_count;
    }
    for (unsigned int corner = 0; corner < 8; ++corner)
    {
      const auto index = data.corners[n.first][corner];
      if (!require(index < data.values.size(), error, "invalid certified corner reference"))
      {
        return false;
      }
      used_values[index] = 1;
      const auto& v = data.values[index];
      const auto coordinate = Training::CornerCoordinate(bounds, corner);
      const double radius = coordinate[1];
      if (!require(radius + v[0] >= data.geometry.inner_radius_cm - Training::EndpointRadialToleranceCm &&
                       radius + v[0] <= data.geometry.outer_radius_cm + Training::EndpointRadialToleranceCm && std::abs(v[1] / radius) <= .5 * TwoPi &&
                       (std::abs(coordinate[2]) != data.geometry.readout_abs_z_cm || v == Value{{0., 0., 0.}}),
                   error, "invalid certified endpoint/readout/angular value"))
      {
        return false;
      }
    }
    return require(Training::UsableCornerValues(data.geometry, bounds, data.maximum_drift_time_ns,
                                                [&](unsigned int corner) -> const Value*
                                                { return &data.values[data.corners[n.first][corner]]; }),
                   error, "unsafe corner angular branch or physical endpoint");
  };
  for (std::uint64_t root = 0; root < roots; ++root)
  {
    if (root < data.first_root || root >= data.end_root)
    {
      if (!require(data.root_lookup[root] == Invalid, error, "out-of-range compact root populated"))
      {
        return false;
      }
      continue;
    }
    const auto before = visited_count;
    if (!visit(data.root_lookup[root], 0, Training::RootBounds(data.geometry, root)))
    {
      return false;
    }
    if (!require(visited_count - before <= data.max_nodes_per_root, error, "source per-root node budget exceeded"))
    {
      return false;
    }
  }
  return require(certified_count == proof.certified_leaves_checked && uncertified_count == proof.uncertified_leaves_checked &&
                     std::all_of(visited.begin(), visited.end(), [](unsigned char v)
                                 { return v == 1; }) &&
                     std::all_of(used_corners.begin(), used_corners.end(), [](unsigned char v)
                                 { return v == 1; }) &&
                     std::all_of(used_values.begin(), used_values.end(), [](unsigned char v)
                                 { return v == 1; }),
                 error, "orphan compact node/corner/value");
}

bool PHGarfieldCompactMap::FindLeaf(const Data& data, unsigned int side, double phi, double r, double z, Lookup& result)
{
  result = {};
  const auto& g = data.geometry;
  if (side > 1 || !std::isfinite(phi) || !std::isfinite(r) || !std::isfinite(z) || r < g.inner_radius_cm || r > g.outer_radius_cm ||
      (side ? z < g.central_membrane_half_width_cm || z > g.readout_abs_z_cm : z < -g.readout_abs_z_cm || z > -g.central_membrane_half_width_cm))
  {
    return false;
  }
  phi = Training::CanonicalPhi(phi);
  const auto root = Training::RootIndex(g, side, phi, r, z);
  if (root >= data.root_lookup.size())
  {
    return false;
  }
  auto id = data.root_lookup[root];
  auto bounds = Training::RootBounds(g, root);
  const std::array<double, 3> coordinate{{phi, r, z}};
  for (unsigned int depth = 0; depth <= data.max_depth; ++depth)
  {
    if (id >= data.nodes.size())
    {
      return false;
    }
    const auto& n = data.nodes[id];
    if (kind(n) == 3)
    {
      result = {id, n.first, depth, leaf_class(n), reason(n), bounds};
      return true;
    }
    const unsigned int axis = kind(n);
    const unsigned int child = coordinate[axis] >= .5 * (bounds.lo[axis] + bounds.hi[axis]);
    id = child ? n.second : n.first;
    bounds = Training::ChildBounds(bounds, axis, child);
  }
  return false;
}
bool PHGarfieldCompactMap::Interpolate(const Data& data, const Lookup& cell, double phi, double r, double z, PHGarfieldForwardMap::NodeValue& result)
{
  result = {};
  if (cell.leaf_class != static_cast<int>(LeafClass::Certified) || cell.corner_record >= data.corners.size())
  {
    return false;
  }
  return Training::InterpolateCorners(cell.bounds, {{Training::CanonicalPhi(phi), r, z}}, [&](unsigned int corner) -> const Value*
                                      {
                                        const auto id = data.corners[cell.corner_record][corner];
                                        return id < data.values.size() ? &data.values[id] : nullptr; }, result);
}

bool PHGarfieldCompactMap::HasUsableCorners(const Data& data, const Lookup& cell)
{
  const bool certified = cell.leaf_class == static_cast<int>(LeafClass::Certified);
  if (cell.cell_id >= data.nodes.size() || cell.corner_record >= data.corners.size() ||
      data.nodes[cell.cell_id].first != cell.corner_record ||
      (!certified && (data.schema_version != SchemaVersion ||
                      (cell.leaf_class != static_cast<int>(LeafClass::Uncertified) &&
                       cell.leaf_class != static_cast<int>(LeafClass::NumericallyUnresolved)) ||
                      cell.fallback_reason == static_cast<int>(FallbackReason::AngularBranch))))
  {
    return false;
  }
  return Training::UsableCornerValues(data.geometry, cell.bounds, data.maximum_drift_time_ns,
                                      [&](unsigned int corner) -> const Value*
                                      {
                                        const auto id = data.corners[cell.corner_record][corner];
                                        return id < data.values.size() ? &data.values[id] : nullptr;
                                      });
}

bool PHGarfieldCompactMap::InterpolateRelaxed(const Data& data, const Lookup& cell, double phi, double r, double z,
                                              PHGarfieldForwardMap::NodeValue& result, bool& certified)
{
  result = {};
  certified = false;
  if (cell.leaf_class == static_cast<int>(LeafClass::Certified))
  {
    certified = Interpolate(data, cell, phi, r, z, result);
    return certified;
  }
  return HasUsableCorners(data, cell) &&
         Training::InterpolateCorners(cell.bounds, {{Training::CanonicalPhi(phi), r, z}}, [&](unsigned int corner) -> const Value*
                                      { return &data.values[data.corners[cell.corner_record][corner]]; }, result);
}

std::string PHGarfieldCompactMap::ContentFingerprint(const Data& data)
{
  PHGarfieldFingerprint::Fnv1a64 hash;
  const auto integer = [&](std::uint64_t value)
  { hash.Integer(value); };
  const auto real = [&](double value)
  { hash.Real(value); };
  const auto text = [&](const std::string& value)
  {
    hash.Integer(value.size());
    hash.Bytes(value);
  };
  const auto* layout = LayoutForSchema(data.schema_version);
  text(layout ? layout : "unsupported_compact_layout");
  integer(data.geometry.phi_nodes);
  integer(data.geometry.radial_nodes);
  integer(data.geometry.longitudinal_nodes);
  real(data.geometry.inner_radius_cm);
  real(data.geometry.outer_radius_cm);
  real(data.geometry.central_membrane_half_width_cm);
  real(data.geometry.readout_abs_z_cm);
  real(data.geometry.nominal_velocity_cm_ns);
  text(data.source.file_sha256);
  text(data.source.configuration_fingerprint_fnv1a64);
  text(data.source.content_fingerprint_fnv1a64);
  text(data.source.canonical_policy);
  text(data.source.converter_sha256);
  integer(data.first_root);
  integer(data.end_root);
  integer(data.max_depth);
  integer(data.max_nodes_per_root);
  integer(data.max_queries_per_root);
  real(data.guard_um);
  real(data.maximum_drift_time_ns);
  integer(data.nodes.size());
  integer(data.corners.size());
  integer(data.values.size());
  integer(data.root_lookup.size());
  for (const auto& n : data.nodes)
  {
    integer(n.first);
    integer(n.second);
    integer(n.tag);
  }
  for (const auto& corners : data.corners)
  {
    for (const auto index : corners)
    {
      integer(index);
    }
  }
  for (const auto& value : data.values)
  {
    for (const auto component : value)
    {
      real(component);
    }
  }
  for (const auto root : data.root_lookup)
  {
    integer(root);
  }
  return hash.Hex();
}
std::string PHGarfieldCompactMap::ProofCanonical(const Data& data)
{
  const auto& p = data.proof;
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << (data.schema_version == LegacySchemaVersion ? LegacyConversionContract : ConversionContract)
      << ";source_sha256=" << data.source.file_sha256 << ";source_config=" << data.source.configuration_fingerprint_fnv1a64
      << ";source_content=" << data.source.content_fingerprint_fnv1a64 << ";converter_sha256=" << data.source.converter_sha256
      << ";compact_content=" << ContentFingerprint(data) << ";source_points=" << p.source_points << ";source_cells=" << p.source_cells
      << ";roots=" << p.roots_checked << ";bounds=" << p.bounds_checked << ";routes=" << p.routes_checked << ";certified=" << p.certified_leaves_checked
      << ";corner_values=" << p.corner_values_checked;
  if (data.schema_version == SchemaVersion)
  {
    out << ";uncertified=" << p.uncertified_leaves_checked;
  }
  return out.str();
}

namespace
{
  using Field = PHGarfieldPrimitiveTree::Field;
  template <class T>
  std::vector<Field> fields(T&);
  template <>
  std::vector<Field> fields(Map::Node& n)
  {
    return {{"first", &n.first, "UInt_t", 1, "first/i"}, {"second", &n.second, "UInt_t", 1, "second/i"}, {"tag", &n.tag, "UInt_t", 1, "tag/i"}};
  }
  template <>
  std::vector<Field> fields(Map::Corners& c)
  {
    return {{"corner_id", c.data(), "UInt_t", 8, "corner_id[8]/i"}};
  }
  template <>
  std::vector<Field> fields(Map::Value& v)
  {
    return {{"value", v.data(), "Double_t", 3, "value[3]/D"}};
  }
  template <>
  std::vector<Field> fields(std::uint32_t& r)
  {
    return {{"cell_id", &r, "UInt_t", 1, "cell_id/i"}};
  }
  template <class T>
  bool write_tree(TDirectory* directory, const char* name, const std::vector<T>& records, std::ostream& error)
  {
    return PHGarfieldPrimitiveTree::Write(directory, name, records, error, fields<T>, "Compact");
  }
  template <class T>
  bool read_tree(TTree* tree, std::vector<T>& records, std::ostream& error)
  {
    return PHGarfieldPrimitiveTree::Read(tree, records, error, fields<T>, "Compact");
  }
  template <class T>
  bool parameter(TDirectory* directory, const char* name, T& value)
  {
    auto* object = dynamic_cast<TParameter<T>*>(directory->Get(name));
    if (!object)
    {
      return false;
    }
    value = object->GetVal();
    return true;
  }
  template <class T>
  bool write_parameter(const char* name, T value)
  {
    TParameter<T> p(name, value);
    return p.Write() > 0;
  }
  bool text(TDirectory* directory, const char* name, std::string& value)
  {
    auto* object = dynamic_cast<TNamed*>(directory->Get(name));
    if (!object)
    {
      return false;
    }
    value = object->GetTitle();
    return true;
  }
}  // namespace

bool PHGarfieldCompactMap::Write(TDirectory* directory, const Data& data, const Limits& limits, std::ostream& error)
{
  if (!directory || !Validate(data, limits, error, false))
  {
    return false;
  }
  auto* file = directory->GetFile();
  if (!require(file && file->IsOpen() && file->IsWritable() && !file->TestBit(TFile::kWriteError), error, "compact output is not writable"))
  {
    return false;
  }
  TDirectory::TContext context(directory);
  const auto& b = data.source;
  const auto& p = data.proof;
  const bool ok = TNamed("layout", LayoutForSchema(data.schema_version)).Write() > 0 &&
                  TNamed("conversion_contract", data.schema_version == LegacySchemaVersion ? LegacyConversionContract : ConversionContract).Write() > 0 &&
                  TNamed("source_file_sha256", b.file_sha256.c_str()).Write() > 0 && TNamed("source_configuration_fingerprint_fnv1a64", b.configuration_fingerprint_fnv1a64.c_str()).Write() > 0 &&
                  TNamed("source_content_fingerprint_fnv1a64", b.content_fingerprint_fnv1a64.c_str()).Write() > 0 && TNamed("source_canonical_policy", b.canonical_policy.c_str()).Write() > 0 &&
                  TNamed("converter_sha256", b.converter_sha256.c_str()).Write() > 0 && TNamed("conversion_proof_sha256", p.sha256.c_str()).Write() > 0 &&
                  TNamed("conversion_proof_canonical", ProofCanonical(data).c_str()).Write() > 0 && TNamed("content_fingerprint_fnv1a64", ContentFingerprint(data).c_str()).Write() > 0 &&
                  write_parameter("schema_version", data.schema_version) && write_parameter<Long64_t>("first_root", data.first_root) && write_parameter<Long64_t>("end_root", data.end_root) &&
                  write_parameter<Long64_t>("node_count", data.nodes.size()) && write_parameter<Long64_t>("corner_record_count", data.corners.size()) && write_parameter<Long64_t>("value_count", data.values.size()) &&
                  write_parameter<Long64_t>("root_lookup_count", data.root_lookup.size()) && write_parameter<Long64_t>("max_depth", data.max_depth) &&
                  write_parameter<Long64_t>("max_nodes_per_root", data.max_nodes_per_root) && write_parameter<Long64_t>("max_queries_per_root", data.max_queries_per_root) &&
                  write_parameter("guard_um", data.guard_um) && write_parameter("maximum_drift_time_ns", data.maximum_drift_time_ns) &&
                  write_parameter<Long64_t>("source_point_count", p.source_points) && write_parameter<Long64_t>("source_cell_count", p.source_cells) &&
                  write_parameter<Long64_t>("proof_roots_checked", p.roots_checked) && write_parameter<Long64_t>("proof_bounds_checked", p.bounds_checked) &&
                  write_parameter<Long64_t>("proof_routes_checked", p.routes_checked) && write_parameter<Long64_t>("proof_certified_leaves_checked", p.certified_leaves_checked) &&
                  write_parameter<Long64_t>("proof_corner_values_checked", p.corner_values_checked);
  if (!require(ok && !file->TestBit(TFile::kWriteError), error, "compact metadata write failed"))
  {
    return false;
  }
  if (data.schema_version == SchemaVersion &&
      !require(write_parameter<Long64_t>("proof_uncertified_leaves_checked", p.uncertified_leaves_checked),
               error, "uncertified corner proof write failed"))
  {
    return false;
  }
  return write_tree(directory, "Nodes", data.nodes, error) && write_tree(directory, "Corners", data.corners, error) &&
         write_tree(directory, "Values", data.values, error) && write_tree(directory, "Roots", data.root_lookup, error);
}

bool PHGarfieldCompactMap::Read(TDirectory* directory, const Geometry& geometry, const ExpectedSource& expected, Data& result, const Limits& limits, std::ostream& error, bool complete)
{
  result = {};
  Data data;
  data.geometry = geometry;
  if (!ValidateLimits(limits, error))
  {
    return false;
  }
  if (!require(directory && Training::RootCount(geometry) && hex(expected.file_sha256, 64) &&
                   hex(expected.configuration_fingerprint_fnv1a64, 16) && hex(expected.content_fingerprint_fnv1a64, 16),
               error, "invalid directory/geometry/external source expectations"))
  {
    return false;
  }
  int version = 0;
  std::string layout;
  std::string contract;
  std::string stored_content;
  std::string stored_proof;
  auto& b = data.source;
  auto& p = data.proof;
  if (!require(parameter(directory, "schema_version", version) && (version == LegacySchemaVersion || version == SchemaVersion) &&
                   text(directory, "layout", layout) && layout == LayoutForSchema(version) &&
                   text(directory, "conversion_contract", contract) && contract == (version == LegacySchemaVersion ? LegacyConversionContract : ConversionContract) && text(directory, "source_file_sha256", b.file_sha256) &&
                   text(directory, "source_configuration_fingerprint_fnv1a64", b.configuration_fingerprint_fnv1a64) && text(directory, "source_content_fingerprint_fnv1a64", b.content_fingerprint_fnv1a64) &&
                   text(directory, "source_canonical_policy", b.canonical_policy) && text(directory, "converter_sha256", b.converter_sha256) &&
                   text(directory, "conversion_proof_sha256", p.sha256) && text(directory, "conversion_proof_canonical", stored_proof) && text(directory, "content_fingerprint_fnv1a64", stored_content) &&
                   binding_valid(b) && b.file_sha256 == expected.file_sha256 && b.configuration_fingerprint_fnv1a64 == expected.configuration_fingerprint_fnv1a64 && b.content_fingerprint_fnv1a64 == expected.content_fingerprint_fnv1a64,
               error, "compact format or external source binding differs"))
  {
    return false;
  }
  data.schema_version = version;
  Long64_t first = 0;
  Long64_t end = 0;
  Long64_t nodes = 0;
  Long64_t corners = 0;
  Long64_t values = 0;
  Long64_t roots = 0;
  Long64_t depth = 0;
  Long64_t per_root = 0;
  Long64_t queries = 0;
  Long64_t source_points = 0;
  Long64_t source_cells = 0;
  Long64_t proof_roots = 0;
  Long64_t proof_bounds = 0;
  Long64_t proof_routes = 0;
  Long64_t proof_certified = 0;
  Long64_t proof_uncertified = 0;
  Long64_t proof_corners = 0;
  if (!require(parameter(directory, "first_root", first) && parameter(directory, "end_root", end) && parameter(directory, "node_count", nodes) &&
                   parameter(directory, "corner_record_count", corners) && parameter(directory, "value_count", values) && parameter(directory, "root_lookup_count", roots) &&
                   parameter(directory, "max_depth", depth) && parameter(directory, "max_nodes_per_root", per_root) && parameter(directory, "max_queries_per_root", queries) &&
                   parameter(directory, "guard_um", data.guard_um) && parameter(directory, "maximum_drift_time_ns", data.maximum_drift_time_ns) &&
                   parameter(directory, "source_point_count", source_points) && parameter(directory, "source_cell_count", source_cells) &&
                   parameter(directory, "proof_roots_checked", proof_roots) && parameter(directory, "proof_bounds_checked", proof_bounds) && parameter(directory, "proof_routes_checked", proof_routes) &&
                   parameter(directory, "proof_certified_leaves_checked", proof_certified) && parameter(directory, "proof_corner_values_checked", proof_corners),
               error, "missing/incorrectly typed compact metadata"))
  {
    return false;
  }
  if (version == SchemaVersion &&
      !require(parameter(directory, "proof_uncertified_leaves_checked", proof_uncertified),
               error, "missing uncertified corner proof"))
  {
    return false;
  }
  std::uint64_t bytes = 0;
  const auto total = Training::RootCount(geometry);
  if (!require(first >= 0 && end >= first && static_cast<std::uint64_t>(end) <= total && (!complete || (!first && static_cast<std::uint64_t>(end) == total)) &&
                   nodes >= 0 && static_cast<std::uint64_t>(nodes) <= limits.max_nodes && corners >= 0 && static_cast<std::uint64_t>(corners) <= limits.max_corners &&
                   values >= 0 && static_cast<std::uint64_t>(values) <= limits.max_values && roots >= 0 && static_cast<std::uint64_t>(roots) == total && total <= limits.max_roots &&
                   depth >= 0 && depth <= 63 && per_root > 0 && per_root <= std::numeric_limits<unsigned int>::max() && queries > 0 && queries <= std::numeric_limits<unsigned int>::max() &&
                   source_points >= values && static_cast<std::uint64_t>(source_points) <= limits.max_source_points && source_cells == nodes && proof_roots == end - first && proof_bounds == nodes && proof_routes == nodes &&
                   proof_certified >= 0 && proof_uncertified >= 0 && proof_certified <= corners && proof_uncertified <= corners - proof_certified &&
                   proof_certified + proof_uncertified == corners && proof_corners == 8 * corners && EstimatedBytes(nodes, corners, values, roots, bytes) && bytes <= limits.max_bytes,
               error, "negative/overflowing/unbounded compact declarations"))
  {
    return false;
  }
  auto* nt = dynamic_cast<TTree*>(directory->Get("Nodes"));
  auto* ct = dynamic_cast<TTree*>(directory->Get("Corners"));
  auto* vt = dynamic_cast<TTree*>(directory->Get("Values"));
  auto* rt = dynamic_cast<TTree*>(directory->Get("Roots"));
  if (!require(nt && ct && vt && rt && nt->GetEntries() == nodes && ct->GetEntries() == corners && vt->GetEntries() == values && rt->GetEntries() == roots, error, "compact tree/declaration row mismatch"))
  {
    return false;
  }
  data.first_root = first;
  data.end_root = end;
  data.max_depth = depth;
  data.max_nodes_per_root = per_root;
  data.max_queries_per_root = queries;
  p.source_points = source_points;
  p.source_cells = source_cells;
  p.roots_checked = proof_roots;
  p.bounds_checked = proof_bounds;
  p.routes_checked = proof_routes;
  p.certified_leaves_checked = proof_certified;
  p.uncertified_leaves_checked = proof_uncertified;
  p.corner_values_checked = proof_corners;
  data.nodes.resize(nodes);
  data.corners.resize(corners);
  data.values.resize(values);
  data.root_lookup.resize(roots);
  if (!read_tree(nt, data.nodes, error) || !read_tree(ct, data.corners, error) || !read_tree(vt, data.values, error) || !read_tree(rt, data.root_lookup, error) ||
      !require(ContentFingerprint(data) == stored_content && ProofCanonical(data) == stored_proof, error, "compact content/conversion proof changed") ||
      !Validate(data, limits, error, complete))
  {
    return false;
  }
  result = std::move(data);
  return true;
}
