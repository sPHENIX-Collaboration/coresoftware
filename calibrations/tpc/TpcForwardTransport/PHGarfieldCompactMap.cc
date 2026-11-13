#include "PHGarfieldCompactMap.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include "PHGarfieldCompactMapInternal.h"
#include "PHGarfieldFingerprint.h"
#include "PHGarfieldPrimitiveTree.h"

namespace
{
  using Map = PHGarfieldCompactMap;
  using Training = PHGarfieldAdaptiveMap;
  using PHGarfieldCompactMapInternal::class_reason;
  using PHGarfieldCompactMapInternal::hex;
  using PHGarfieldCompactMapInternal::kind;
  using PHGarfieldCompactMapInternal::leaf_class;
  using PHGarfieldCompactMapInternal::reason;
  using PHGarfieldCompactMapInternal::require;
}  // namespace

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
  return require(limits.max_nodes && limits.max_nodes <= HardMaximumNodes &&
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
  if (!require(roots && roots <= limits.max_roots && hex(data.configuration_fingerprint_fnv1a64, 16) &&
                   !data.canonical_policy.empty() && data.canonical_policy.size() <= 65536 &&
                   data.first_root <= data.end_root && data.end_root <= roots && (!complete || (!data.first_root && data.end_root == roots)) &&
                   data.root_lookup.size() == roots && data.nodes.size() <= limits.max_nodes && data.corners.size() <= limits.max_corners && data.values.size() <= limits.max_values &&
                   data.nodes.capacity() <= limits.max_nodes && data.corners.capacity() <= limits.max_corners && data.values.capacity() <= limits.max_values && data.root_lookup.capacity() <= limits.max_roots &&
                   EstimatedBytes(data.nodes.capacity(), data.corners.capacity(), data.values.capacity(), data.root_lookup.capacity(), bytes) && bytes <= limits.max_bytes &&
                   data.max_depth <= 63 && data.max_nodes_per_root && data.max_queries_per_root && std::isfinite(data.guard_um) && data.guard_um > 0. &&
                   std::isfinite(data.maximum_drift_time_ns) && data.maximum_drift_time_ns >= 0.,
               error, "invalid geometry/configuration/count/capacity/budget"))
  {
    return false;
  }
  std::vector<unsigned char> visited(data.nodes.size(), 0);
  std::vector<unsigned char> used_corners(data.corners.size(), 0);
  std::vector<unsigned char> used_values(data.values.size(), 0);
  std::uint64_t visited_count = 0;
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
        !require((c == static_cast<int>(LeafClass::Uncertified) || c == static_cast<int>(LeafClass::NumericallyUnresolved)) &&
                     r != static_cast<int>(FallbackReason::AngularBranch),
                 error, "ineligible noncertified leaf has interpolation payload"))
    {
      return false;
    }
    if (!require(n.first < data.corners.size() && !used_corners[n.first], error, "invalid/duplicate corner record"))
    {
      return false;
    }
    used_corners[n.first] = 1;
    for (unsigned int corner = 0; corner < 8; ++corner)
    {
      const auto index = data.corners[n.first][corner];
      if (!require(index < data.values.size(), error, "invalid certified corner reference"))
      {
        return false;
      }
      used_values[index] = 1;
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
  return require(std::all_of(visited.begin(), visited.end(), [](unsigned char v)
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
      (!certified && ((cell.leaf_class != static_cast<int>(LeafClass::Uncertified) &&
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
  integer(SchemaVersion);
  text(Layout);
  integer(data.geometry.phi_nodes);
  integer(data.geometry.radial_nodes);
  integer(data.geometry.longitudinal_nodes);
  real(data.geometry.inner_radius_cm);
  real(data.geometry.outer_radius_cm);
  real(data.geometry.central_membrane_half_width_cm);
  real(data.geometry.readout_abs_z_cm);
  real(data.geometry.nominal_velocity_cm_ns);
  text(data.configuration_fingerprint_fnv1a64);
  text(data.canonical_policy);
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
  const bool ok = TNamed("layout", Layout).Write() > 0 &&
                  TNamed("configuration_fingerprint_fnv1a64", data.configuration_fingerprint_fnv1a64.c_str()).Write() > 0 &&
                  TNamed("canonical_policy", data.canonical_policy.c_str()).Write() > 0 &&
                  TNamed("content_fingerprint_fnv1a64", ContentFingerprint(data).c_str()).Write() > 0 &&
                  write_parameter("schema_version", SchemaVersion) && write_parameter<Long64_t>("first_root", data.first_root) && write_parameter<Long64_t>("end_root", data.end_root) &&
                  write_parameter<Long64_t>("node_count", data.nodes.size()) && write_parameter<Long64_t>("corner_record_count", data.corners.size()) && write_parameter<Long64_t>("value_count", data.values.size()) &&
                  write_parameter<Long64_t>("root_lookup_count", data.root_lookup.size()) && write_parameter<Long64_t>("max_depth", data.max_depth) &&
                  write_parameter<Long64_t>("max_nodes_per_root", data.max_nodes_per_root) && write_parameter<Long64_t>("max_queries_per_root", data.max_queries_per_root) &&
                  write_parameter("guard_um", data.guard_um) && write_parameter("maximum_drift_time_ns", data.maximum_drift_time_ns);
  if (!require(ok && !file->TestBit(TFile::kWriteError), error, "compact metadata write failed"))
  {
    return false;
  }
  return write_tree(directory, "Nodes", data.nodes, error) && write_tree(directory, "Corners", data.corners, error) &&
         write_tree(directory, "Values", data.values, error) && write_tree(directory, "Roots", data.root_lookup, error);
}

bool PHGarfieldCompactMap::Read(TDirectory* directory, const Geometry& geometry, Data& result, const Limits& limits, std::ostream& error, bool complete)
{
  result = {};
  Data data;
  data.geometry = geometry;
  const auto total = Training::RootCount(geometry);
  if (!ValidateLimits(limits, error) || !require(directory && total && total <= limits.max_roots, error, "invalid directory/geometry/root budget"))
  {
    return false;
  }
  int version = 0;
  std::string layout;
  std::string stored_content;
  if (!require(parameter(directory, "schema_version", version) && version == SchemaVersion &&
                   text(directory, "layout", layout) && layout == Layout &&
                   text(directory, "configuration_fingerprint_fnv1a64", data.configuration_fingerprint_fnv1a64) &&
                   text(directory, "canonical_policy", data.canonical_policy) && text(directory, "content_fingerprint_fnv1a64", stored_content) &&
                   hex(data.configuration_fingerprint_fnv1a64, 16) && hex(stored_content, 16) &&
                   !data.canonical_policy.empty() && data.canonical_policy.size() <= 65536,
               error, "unsupported format or missing configuration/content fingerprint"))
  {
    return false;
  }
  Long64_t first = 0;
  Long64_t end = 0;
  Long64_t nodes = 0;
  Long64_t corners = 0;
  Long64_t values = 0;
  Long64_t roots = 0;
  Long64_t depth = 0;
  Long64_t per_root = 0;
  Long64_t queries = 0;
  if (!require(parameter(directory, "first_root", first) && parameter(directory, "end_root", end) && parameter(directory, "node_count", nodes) &&
                   parameter(directory, "corner_record_count", corners) && parameter(directory, "value_count", values) && parameter(directory, "root_lookup_count", roots) &&
                   parameter(directory, "max_depth", depth) && parameter(directory, "max_nodes_per_root", per_root) && parameter(directory, "max_queries_per_root", queries) &&
                   parameter(directory, "guard_um", data.guard_um) && parameter(directory, "maximum_drift_time_ns", data.maximum_drift_time_ns),
               error, "missing/incorrectly typed compact metadata"))
  {
    return false;
  }
  std::uint64_t bytes = 0;
  if (!require(first >= 0 && end >= first && static_cast<std::uint64_t>(end) <= total && (!complete || (!first && static_cast<std::uint64_t>(end) == total)) &&
                   nodes >= 0 && static_cast<std::uint64_t>(nodes) <= limits.max_nodes && corners >= 0 && static_cast<std::uint64_t>(corners) <= limits.max_corners &&
                   values >= 0 && static_cast<std::uint64_t>(values) <= limits.max_values && roots >= 0 && static_cast<std::uint64_t>(roots) == total &&
                   depth >= 0 && depth <= 63 && per_root > 0 && per_root <= std::numeric_limits<unsigned int>::max() && queries > 0 && queries <= std::numeric_limits<unsigned int>::max() &&
                   std::isfinite(data.guard_um) && data.guard_um > 0. && std::isfinite(data.maximum_drift_time_ns) && data.maximum_drift_time_ns >= 0. &&
                   EstimatedBytes(nodes, corners, values, roots, bytes) && bytes <= limits.max_bytes,
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
  data.nodes.resize(nodes);
  data.corners.resize(corners);
  data.values.resize(values);
  data.root_lookup.resize(roots);
  if (!read_tree(nt, data.nodes, error) || !read_tree(ct, data.corners, error) || !read_tree(vt, data.values, error) || !read_tree(rt, data.root_lookup, error) ||
      !require(ContentFingerprint(data) == stored_content, error, "compact content changed") || !Validate(data, limits, error, complete))
  {
    return false;
  }
  result = std::move(data);
  return true;
}
