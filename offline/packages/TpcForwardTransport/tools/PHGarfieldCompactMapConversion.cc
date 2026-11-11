#include "PHGarfieldCompactMap.h"
#include "PHGarfieldCompactMapInternal.h"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <ostream>
#include <utility>
#include <vector>

namespace
{
  using Training = PHGarfieldAdaptiveMap;
  using PHGarfieldCompactMapInternal::binding_valid;
  using PHGarfieldCompactMapInternal::bits;
  using PHGarfieldCompactMapInternal::bounds_equal;
  using PHGarfieldCompactMapInternal::geometry_equal;
  using PHGarfieldCompactMapInternal::kind;
  using PHGarfieldCompactMapInternal::leaf_class;
  using PHGarfieldCompactMapInternal::reason;
  using PHGarfieldCompactMapInternal::require;
}  // namespace

bool PHGarfieldCompactMap::VerifyConversion(const Training::Data& source, const Data& data, std::ostream& error)
{
  if (!require(data.proof.source_points == source.points.size() && data.proof.source_cells == source.cells.size() &&
                   data.proof.roots_checked == source.roots.size() && data.proof.bounds_checked == source.cells.size() && data.proof.routes_checked == source.cells.size() &&
                   geometry_equal(source.geometry, data.geometry) && source.cells.size() == data.nodes.size() && source.first_root == data.first_root && source.end_root == data.end_root &&
                   source.roots.size() == data.end_root - data.first_root && data.root_lookup.size() == Training::RootCount(data.geometry) &&
                   source.policy.max_depth == data.max_depth && source.policy.max_nodes_per_root == data.max_nodes_per_root && source.policy.max_queries_per_root == data.max_queries_per_root && bits(source.policy.guard_um, data.guard_um) &&
                   bits(source.maximum_drift_time_ns, data.maximum_drift_time_ns) && data.source.content_fingerprint_fnv1a64 == Training::ContentFingerprint(source) &&
                   data.source.canonical_policy == Training::BuildCanonicalPolicy(source),
               error, "source geometry/policy/content mismatch"))
  {
    return false;
  }
  std::vector<std::uint32_t> point_mapping(source.points.size(), Invalid);
  std::vector<std::uint32_t> back_mapping(data.values.size(), Invalid);
  std::vector<unsigned char> visited(data.nodes.size(), 0);
  std::uint64_t certified_checked = 0;
  std::uint64_t uncertified_checked = 0;
  std::function<bool(std::uint32_t, unsigned int, const Bounds&)> check;
  check = [&](std::uint32_t id, unsigned int depth, const Bounds& bounds)
  {
    if (!require(id < source.cells.size() && id < data.nodes.size() && !visited[id], error, "conversion node index out of range/overlapping"))
    {
      return false;
    }
    visited[id] = 1;
    const auto& original = source.cells.at(id);
    const auto& n = data.nodes.at(id);
    if (!require(original.id == id && original.depth == depth && bounds_equal(Bounds{original.lo, original.hi}, bounds), error, "conversion changed exact bounds/source ID/depth"))
    {
      return false;
    }
    if (original.split_axis >= 0)
    {
      return require(n.tag == static_cast<unsigned int>(original.split_axis) && n.first == original.child_id[0] && n.second == original.child_id[1], error, "conversion changed child/axis routing") &&
             check(n.first, depth + 1, Training::ChildBounds(bounds, original.split_axis, 0)) && check(n.second, depth + 1, Training::ChildBounds(bounds, original.split_axis, 1));
    }
    if (!require(kind(n) == 3 && leaf_class(n) == original.leaf_class && reason(n) == original.fallback_reason, error, "conversion changed leaf class/reason"))
    {
      return false;
    }
    const bool retain = original.leaf_class == static_cast<int>(LeafClass::Certified) ||
                        (data.schema_version == SchemaVersion && Training::HasUsableCorners(source, original));
    if (!retain)
    {
      return require(n.first == Invalid && n.second == 0, error, "direct leaf acquired corner payload");
    }
    if (!require(n.first < data.corners.size(), error, "conversion missing corner record"))
    {
      return false;
    }
    if (original.leaf_class == static_cast<int>(LeafClass::Certified))
    {
      ++certified_checked;
    }
    else
    {
      ++uncertified_checked;
    }
    for (unsigned int corner = 0; corner < 8; ++corner)
    {
      const auto old_id = original.corner_id[corner];
      const auto new_id = data.corners[n.first][corner];
      if (!require(old_id >= 0 && static_cast<std::uint64_t>(old_id) < source.points.size() && new_id < data.values.size(), error, "conversion corner out of bounds"))
      {
        return false;
      }
      if (point_mapping[old_id] == Invalid)
      {
        point_mapping[old_id] = new_id;
      }
      if (back_mapping[new_id] == Invalid)
      {
        back_mapping[new_id] = old_id;
      }
      if (!require(point_mapping[old_id] == new_id && back_mapping[new_id] == old_id, error, "conversion changed point identity deduplication"))
      {
        return false;
      }
      for (unsigned int component = 0; component < 3; ++component)
      {
        if (!require(bits(source.points[old_id].value[component], data.values[new_id][component]), error, "conversion changed corner double bits"))
        {
          return false;
        }
      }
    }
    return true;
  };
  for (const auto& root : source.roots)
  {
    if (!require(root.root_id < data.root_lookup.size() && data.root_lookup[root.root_id] == root.cell_id, error, "conversion changed root lookup") ||
        !check(root.cell_id, 0, Training::RootBounds(data.geometry, root.root_id)))
    {
      return false;
    }
  }
  return require(certified_checked == data.proof.certified_leaves_checked &&
                     uncertified_checked == data.proof.uncertified_leaves_checked &&
                     data.proof.corner_values_checked == 8ULL * (certified_checked + uncertified_checked) &&
                     std::all_of(visited.begin(), visited.end(), [](unsigned char value)
                                 { return value == 1; }) &&
                     std::all_of(back_mapping.begin(), back_mapping.end(), [](std::uint32_t value)
                                 { return value != Invalid; }),
                 error, "conversion proof omitted nodes/values");
}

bool PHGarfieldCompactMap::Convert(Training::Data& source, const SourceBinding& binding, Data& result, const Limits& limits,
                                   std::ostream& error, bool complete, bool retain_uncertified_corners)
{
  result = {};
  if (!ValidateLimits(limits, error))
  {
    return false;
  }
  if (!require(binding_valid(binding) && source.points.size() <= limits.max_source_points && source.cells.size() <= limits.max_nodes &&
                   binding.content_fingerprint_fnv1a64 == Training::ContentFingerprint(source) && binding.canonical_policy == Training::BuildCanonicalPolicy(source),
               error, "invalid source binding/counts") ||
      !Training::Validate(source, error, complete))
  {
    return false;
  }
  Data data;
  data.schema_version = retain_uncertified_corners ? SchemaVersion : LegacySchemaVersion;
  data.geometry = source.geometry;
  data.source = binding;
  data.first_root = source.first_root;
  data.end_root = source.end_root;
  data.max_depth = source.policy.max_depth;
  data.max_nodes_per_root = source.policy.max_nodes_per_root;
  data.max_queries_per_root = source.policy.max_queries_per_root;
  data.guard_um = source.policy.guard_um;
  data.maximum_drift_time_ns = source.maximum_drift_time_ns;
  const auto roots = Training::RootCount(data.geometry);
  if (!require(roots <= limits.max_roots, error, "too many compact roots"))
  {
    return false;
  }
  const auto retain = [&](const Training::Cell& cell)
  {
    return cell.split_axis == -1 && (cell.leaf_class == static_cast<int>(LeafClass::Certified) ||
                                     (retain_uncertified_corners && Training::HasUsableCorners(source, cell)));
  };
  std::uint64_t certified = 0;
  std::uint64_t retained = 0;
  for (const auto& c : source.cells)
  {
    if (c.split_axis == -1 && c.leaf_class == static_cast<int>(LeafClass::Certified))
    {
      ++certified;
    }
    retained += retain(c);
  }
  if (!require(retained <= limits.max_corners, error, "too many retained leaf corner records"))
  {
    return false;
  }
  // The mark array is conversion-only; runtime never stores training IDs.
  std::vector<std::uint32_t> mapping(source.points.size(), Invalid);
  std::uint64_t value_count = 0;
  for (const auto& c : source.cells)
  {
    if (retain(c))
    {
      for (const auto id : c.corner_id)
      {
        if (mapping[id] == Invalid)
        {
          mapping[id] = 0;
          ++value_count;
        }
      }
    }
  }
  std::uint64_t bytes = 0;
  if (!require(value_count <= limits.max_values && EstimatedBytes(source.cells.size(), retained, value_count, roots, bytes) && bytes <= limits.max_bytes,
               error, "compact result exceeds declared allocation budget"))
  {
    return false;
  }
  data.nodes.resize(source.cells.size());
  data.corners.resize(retained);
  data.values.resize(value_count);
  data.root_lookup.assign(roots, Invalid);
  std::fill(mapping.begin(), mapping.end(), Invalid);
  std::uint32_t next_value = 0;
  std::uint32_t next_corner = 0;
  for (const auto& c : source.cells)
  {
    auto& n = data.nodes[c.id];
    if (c.split_axis >= 0)
    {
      n = {static_cast<std::uint32_t>(c.child_id[0]), static_cast<std::uint32_t>(c.child_id[1]), static_cast<std::uint32_t>(c.split_axis)};
      continue;
    }
    n.tag = 3U | (static_cast<unsigned int>(c.leaf_class) << 2U) | (static_cast<unsigned int>(c.fallback_reason) << 5U);
    n.second = 0;
    n.first = Invalid;
    if (!retain(c))
    {
      continue;
    }
    n.first = next_corner++;
    auto& corners = data.corners[n.first];
    for (unsigned int corner = 0; corner < 8; ++corner)
    {
      const auto id = c.corner_id[corner];
      if (mapping[id] == Invalid)
      {
        mapping[id] = next_value;
        data.values[next_value++] = source.points[id].value;
      }
      corners[corner] = mapping[id];
    }
  }
  for (const auto& root : source.roots)
  {
    data.root_lookup[root.root_id] = root.cell_id;
  }
  data.proof = {source.points.size(), source.cells.size(), source.roots.size(), source.cells.size(), source.cells.size(), certified, 8ULL * retained, retained - certified, ""};
  data.proof.sha256 = Sha256(ProofCanonical(data));
  if (!Validate(data, limits, error, complete) || !VerifyConversion(source, data, error))
  {
    return false;
  }
  result = std::move(data);
  return true;
}
