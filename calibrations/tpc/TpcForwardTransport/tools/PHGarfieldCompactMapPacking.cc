#include "PHGarfieldCompactMap.h"
#include "PHGarfieldCompactMapInternal.h"

#include <algorithm>
#include <cstdint>
#include <ostream>
#include <utility>
#include <vector>

namespace
{
  using Training = PHGarfieldAdaptiveMap;
  using PHGarfieldCompactMapInternal::geometry_equal;
  using PHGarfieldCompactMapInternal::hex;
  using PHGarfieldCompactMapInternal::kind;
  using PHGarfieldCompactMapInternal::require;
}  // namespace

bool PHGarfieldCompactMap::Pack(const Training::Data& source, const std::string& configuration_fingerprint, Data& result, const Limits& limits,
                                std::ostream& error, bool complete)
{
  result = {};
  if (!ValidateLimits(limits, error))
  {
    return false;
  }
  if (!require(hex(configuration_fingerprint, 16) && source.cells.size() <= limits.max_nodes,
               error, "invalid configuration fingerprint or node budget") ||
      !Training::Validate(source, error, complete))
  {
    return false;
  }
  Data data;
  data.geometry = source.geometry;
  data.configuration_fingerprint_fnv1a64 = configuration_fingerprint;
  data.canonical_policy = Training::BuildCanonicalPolicy(source);
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
                                     Training::HasUsableCorners(source, cell));
  };
  std::uint64_t retained = 0;
  for (const auto& c : source.cells)
  {
    retained += retain(c);
  }
  if (!require(retained <= limits.max_corners, error, "too many retained leaf corner records"))
  {
    return false;
  }
  // The mark array is construction-only; runtime stores only retained values.
  if (!require(source.points.size() <= (limits.max_bytes - 1024) / sizeof(std::uint32_t), error, "packing workspace exceeds allocation budget"))
  {
    return false;
  }
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
  if (!Validate(data, limits, error, complete))
  {
    return false;
  }
  result = std::move(data);
  return true;
}

bool PHGarfieldCompactMap::Append(Data& result, const Data& shard, const Limits& limits, std::ostream& error)
{
  if (!require(&result != &shard, error, "cannot append a table to itself") || !Validate(shard, limits, error, false))
  {
    return false;
  }
  if (result.root_lookup.empty())
  {
    if (!require(result.nodes.empty() && result.corners.empty() && result.values.empty(), error, "incomplete merge accumulator"))
    {
      return false;
    }
    result.geometry = shard.geometry;
    result.configuration_fingerprint_fnv1a64 = shard.configuration_fingerprint_fnv1a64;
    result.canonical_policy = shard.canonical_policy;
    result.first_root = result.end_root = shard.first_root;
    result.max_depth = shard.max_depth;
    result.max_nodes_per_root = shard.max_nodes_per_root;
    result.max_queries_per_root = shard.max_queries_per_root;
    result.guard_um = shard.guard_um;
    result.maximum_drift_time_ns = 0.;
    result.root_lookup.assign(shard.root_lookup.size(), Invalid);
  }
  const auto node_offset = result.nodes.size();
  const auto corner_offset = result.corners.size();
  const auto value_offset = result.values.size();
  std::uint64_t bytes = 0;
  if (!require(geometry_equal(result.geometry, shard.geometry) && result.configuration_fingerprint_fnv1a64 == shard.configuration_fingerprint_fnv1a64 &&
                   result.canonical_policy == shard.canonical_policy && result.end_root == shard.first_root &&
                   result.max_depth == shard.max_depth && result.max_nodes_per_root == shard.max_nodes_per_root && result.max_queries_per_root == shard.max_queries_per_root &&
                   PHGarfieldCompactMapInternal::bits(result.guard_um, shard.guard_um) && result.root_lookup.size() == shard.root_lookup.size() &&
                   node_offset <= limits.max_nodes && shard.nodes.size() <= limits.max_nodes - node_offset &&
                   corner_offset <= limits.max_corners && shard.corners.size() <= limits.max_corners - corner_offset &&
                   value_offset <= limits.max_values && shard.values.size() <= limits.max_values - value_offset &&
                   result.nodes.capacity() <= limits.max_nodes && result.corners.capacity() <= limits.max_corners && result.values.capacity() <= limits.max_values && result.root_lookup.capacity() <= limits.max_roots &&
                   EstimatedBytes(node_offset + shard.nodes.size(), corner_offset + shard.corners.size(), value_offset + shard.values.size(), result.root_lookup.size(), bytes) && bytes <= limits.max_bytes,
               error, "incompatible/noncontiguous shards or merged allocation budget exceeded"))
  {
    return false;
  }
  // Geometric capacity growth avoids copying the accumulated table for every
  // shard. Fall back to exact capacities when spare capacity exceeds the budget.
  const auto capacity = [](std::uint64_t current, std::uint64_t needed, std::uint64_t limit)
  { return current >= needed ? current : std::min(limit, std::max(needed, current * 2)); };
  auto nodes_capacity = capacity(result.nodes.capacity(), node_offset + shard.nodes.size(), limits.max_nodes);
  auto corners_capacity = capacity(result.corners.capacity(), corner_offset + shard.corners.size(), limits.max_corners);
  auto values_capacity = capacity(result.values.capacity(), value_offset + shard.values.size(), limits.max_values);
  if (!EstimatedBytes(nodes_capacity, corners_capacity, values_capacity, result.root_lookup.capacity(), bytes) || bytes > limits.max_bytes)
  {
    nodes_capacity = std::max<std::uint64_t>(result.nodes.capacity(), node_offset + shard.nodes.size());
    corners_capacity = std::max<std::uint64_t>(result.corners.capacity(), corner_offset + shard.corners.size());
    values_capacity = std::max<std::uint64_t>(result.values.capacity(), value_offset + shard.values.size());
  }
  if (!require(EstimatedBytes(nodes_capacity, corners_capacity, values_capacity, result.root_lookup.capacity(), bytes) && bytes <= limits.max_bytes,
               error, "merged capacities exceed allocation budget"))
  {
    return false;
  }
  result.nodes.reserve(nodes_capacity);
  result.corners.reserve(corners_capacity);
  result.values.reserve(values_capacity);
  for (auto node : shard.nodes)
  {
    if (kind(node) < 3U)
    {
      node.first += node_offset;
      node.second += node_offset;
    }
    else if (node.first != Invalid)
    {
      node.first += corner_offset;
    }
    result.nodes.push_back(node);
  }
  for (auto corners : shard.corners)
  {
    for (auto& index : corners)
    {
      index += value_offset;
    }
    result.corners.push_back(corners);
  }
  result.values.insert(result.values.end(), shard.values.begin(), shard.values.end());
  for (auto root = shard.first_root; root < shard.end_root; ++root)
  {
    result.root_lookup[root] = shard.root_lookup[root] + node_offset;
  }
  result.end_root = shard.end_root;
  result.maximum_drift_time_ns = std::max(result.maximum_drift_time_ns, shard.maximum_drift_time_ns);
  return true;
}
