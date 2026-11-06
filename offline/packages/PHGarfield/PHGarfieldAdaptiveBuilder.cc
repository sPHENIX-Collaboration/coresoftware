#include "PHGarfieldAdaptiveMap.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <ostream>
#include <queue>
#include <set>

namespace
{
  using Map = PHGarfieldAdaptiveMap;
  constexpr double Pi = 3.1415926535897932384626433832795;

  std::array<double, 2> error_at(const std::array<double, 3>& prediction,
                                 const Map::Point& actual, double velocity)
  {
    const double r = actual.coordinate[1];
    const double pa = r + prediction[0];
    const double aa = prediction[1] / r;
    const double pb = r + actual.value[0];
    const double ab = actual.value[1] / r;
    return {{1.e4 * std::hypot(pa * std::cos(aa) - pb * std::cos(ab), pa * std::sin(aa) - pb * std::sin(ab)),
             1.e4 * velocity * std::abs(prediction[2] - actual.value[2])}};
  }

  std::array<double, 3> linear(const Map::Data& data, const Map::Cell& cell,
                               const std::array<double, 3>& coordinate)
  {
    std::array<double, 3> result{};
    std::array<double, 3> fraction{};
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
      fraction[axis] = (coordinate[axis] - cell.lo[axis]) / (cell.hi[axis] - cell.lo[axis]);
    }
    for (unsigned int corner = 0; corner < 8; ++corner)
    {
      const double weight = ((corner & 4) ? fraction[0] : 1. - fraction[0]) * ((corner & 2) ? fraction[1] : 1. - fraction[1]) * ((corner & 1) ? fraction[2] : 1. - fraction[2]);
      for (unsigned int component = 0; component < 3; ++component)
      {
        result[component] += weight * data.points[cell.corner_id[corner]].value[component];
      }
    }
    return result;
  }
}  // namespace

bool PHGarfieldAdaptiveMap::Build(const Geometry& geometry, const Policy& policy,
                                  std::uint64_t first, std::uint64_t end, const Sampler& sampler, Data& output,
                                  BuildSummary& summary, std::ostream& log, std::ostream& error)
{
  output = {};
  summary = {};
  Data data;
  data.geometry = geometry;
  data.policy = policy;
  data.first_root = first;
  data.end_root = first;
  // Validate policy/global limits before even allocating the coarse index.
  if (!sampler || !Validate(data, error, false) || first >= end || end > RootCount(geometry))
  {
    error << "Adaptive builder: invalid configuration/root range\n";
    return false;
  }
  data.root_lookup.clear();
  data.end_root = end;
  std::map<std::array<double, 3>, std::uint64_t> cache;
  auto memory_available = [&](std::uint64_t points, std::uint64_t cells)
  {
    std::uint64_t bytes = 0;
    const bool ok = points <= policy.max_total_points && cells <= policy.max_total_cells &&
                    EstimatedBytes(points, cells, RootCount(geometry), bytes) && bytes <= policy.max_estimated_bytes;
    if (!ok)
    {
      error << "Adaptive builder: explicit global allocation cap exceeded; no partial LUT may be published\n";
    }
    return ok;
  };
  for (std::uint64_t root = first; root < end; ++root)
  {
    std::set<std::array<double, 3>> root_samples;
    unsigned int root_nodes = 1;
    bool fatal = false;
    // 0=ready,1=per-root sampling budget,2=fatal model/allocation failure.
    auto sample = [&](std::array<double, 3> coordinate, std::int64_t& id)
    {
      coordinate[0] = CanonicalPhi(coordinate[0]);
      if (!root_samples.contains(coordinate))
      {
        if (root_samples.size() >= policy.max_queries_per_root)
        {
          return 1;
        }
        root_samples.insert(coordinate);
      }
      const auto existing = cache.find(coordinate);
      if (existing != cache.end())
      {
        id = existing->second;
        ++summary.reused_points;
        return 0;
      }
      if (!memory_available(data.points.size() + 1, data.cells.size()))
      {
        return 2;
      }
      Point point;
      point.coordinate = coordinate;
      if (!sampler(coordinate, point, error))
      {
        return 2;
      }
      ++summary.actual_sampler_queries;
      point.coordinate = coordinate;
      point.id = data.points.size();
      if (point.reaches)
      {
        data.maximum_drift_time_ns = std::max(data.maximum_drift_time_ns, point.value[2]);
      }
      if (!ValidatePoint(point, geometry, data.maximum_drift_time_ns, error))
      {
        error << "Adaptive builder: invalid sample in root=" << root << " after queries=" << summary.actual_sampler_queries << '\n';
        return 2;
      }
      id = point.id;
      cache.emplace(coordinate, point.id);
      data.points.push_back(point);
      return 0;
    };
    if (!memory_available(data.points.size(), data.cells.size() + 1))
    {
      return false;
    }
    const auto root_cell = data.cells.size();
    Cell initial;
    initial.id = root_cell;
    initial.root_id = root;
    const auto root_bounds = RootBounds(geometry, root);
    initial.lo = root_bounds.lo;
    initial.hi = root_bounds.hi;
    data.cells.push_back(initial);
    data.roots.push_back({root, root_cell});
    struct Pending
    {
      double volume;
      std::uint64_t id;
      bool operator<(const Pending& other) const
      {
        return volume < other.volume || (volume == other.volume && id > other.id);
      }
    };
    std::priority_queue<Pending> pending;
    // Evaluate before scheduling: certified/physical-loss/angular-ambiguous
    // leaves never compete for a refinement quota. The queue favors broad
    // unresolved gas volume; deterministic IDs break equal-volume ties.
    auto evaluate = [&](std::uint64_t id)
    {
      if (fatal)
      {
        return;
      }
      const Bounds bounds{data.cells[id].lo, data.cells[id].hi};
      std::array<std::int64_t, 27> stencil{};
      stencil.fill(-1);
      unsigned int reaching = 0;
      double min_angle = std::numeric_limits<double>::infinity();
      double max_angle = -min_angle;
      // Canonical cache points at phi2pi are phi0, while interpolation fractions
      // below use the declared unwrapped cell stencil coordinates.
      for (unsigned int iphi = 0; iphi < 3; ++iphi)
      {
        for (unsigned int ir = 0; ir < 3; ++ir)
        {
          for (unsigned int iz = 0; iz < 3; ++iz)
          {
            const std::array<unsigned int, 3> index{{iphi, ir, iz}};
            std::array<double, 3> coordinate{};
            for (unsigned int axis = 0; axis < 3; ++axis)
            {
              coordinate[axis] = index[axis] == 0 ? bounds.lo[axis] : (index[axis] == 2 ? bounds.hi[axis] : .5 * (bounds.lo[axis] + bounds.hi[axis]));
            }
            auto& point_id = stencil[(iphi * 3 + ir) * 3 + iz];
            const int status = sample(coordinate, point_id);
            if (status)
            {
              if (status == 2)
              {
                fatal = true;
              }
              data.cells[id].leaf_class = static_cast<int>(LeafClass::Uncertified);
              data.cells[id].fallback_reason = static_cast<int>(FallbackReason::QueryBudget);
              return;
            }
            const auto& point = data.points[point_id];
            reaching += point.reaches;
            const bool corner = iphi != 1 && ir != 1 && iz != 1;
            if (corner)
            {
              data.cells[id].corner_id[((iphi == 2) << 2) | ((ir == 2) << 1) | (iz == 2)] = point_id;
            }
            else
            {
              ++data.cells[id].probe_count;
              data.cells[id].probe_reaching += point.reaches;
            }
            if (point.reaches)
            {
              const double angle = point.value[1] / point.coordinate[1];
              min_angle = std::min(min_angle, angle);
              max_angle = std::max(max_angle, angle);
            }
          }
        }
      }
      auto& c = data.cells[id];
      if (!reaching)
      {
        c.leaf_class = static_cast<int>(LeafClass::AllLost);
        c.fallback_reason = static_cast<int>(FallbackReason::SampledAllLost);
        return;
      }
      if (reaching != 27)
      {
        c.leaf_class = static_cast<int>(LeafClass::Mixed);
        c.fallback_reason = static_cast<int>(FallbackReason::PhysicalMixed);
        for (unsigned int axis = 0; axis < 3; ++axis)
        {
          for (unsigned int other0 = 0; other0 < 3; ++other0)
          {
            for (unsigned int other1 = 0; other1 < 3; ++other1)
            {
              std::array<unsigned int, 3> index{};
              index[(axis + 1) % 3] = other0;
              index[(axis + 2) % 3] = other1;
              std::array<int, 3> reach{};
              for (unsigned int position = 0; position < 3; ++position)
              {
                index[axis] = position;
                reach[position] = data.points[stencil[(index[0] * 3 + index[1]) * 3 + index[2]]].reaches;
              }
              c.directional_reach_changes[axis] += (reach[0] != reach[1]) + (reach[1] != reach[2]);
            }
          }
        }
      }
      else
      {
        c.angular_valid = (max_angle - min_angle <= Pi);
        c.max_transverse_error_um = 0.;
        c.max_time_error_um = 0.;
        for (unsigned int iphi = 0; iphi < 3; ++iphi)
        {
          for (unsigned int ir = 0; ir < 3; ++ir)
          {
            for (unsigned int iz = 0; iz < 3; ++iz)
            {
              const std::array<unsigned int, 3> index{{iphi, ir, iz}};
              std::array<double, 3> coordinate{};
              for (unsigned int axis = 0; axis < 3; ++axis)
              {
                coordinate[axis] = index[axis] == 0 ? bounds.lo[axis] : (index[axis] == 2 ? bounds.hi[axis] : .5 * (bounds.lo[axis] + bounds.hi[axis]));
              }
              const auto difference = error_at(linear(data, c, coordinate), data.points[stencil[(iphi * 3 + ir) * 3 + iz]], geometry.nominal_velocity_cm_ns);
              c.max_transverse_error_um = std::max(c.max_transverse_error_um, difference[0]);
              c.max_time_error_um = std::max(c.max_time_error_um, difference[1]);
            }
          }
        }
        // Each score is measured on nine actual endpoint-midpoint lines parallel
        // to that axis, using all three positions along both other axes. No axis
        // priority or source-map curvature surrogate selects the split.
        for (unsigned int axis = 0; axis < 3; ++axis)
        {
          c.directional_error_um[axis] = 0.;
          for (unsigned int other0 = 0; other0 < 3; ++other0)
          {
            for (unsigned int other1 = 0; other1 < 3; ++other1)
            {
              std::array<unsigned int, 3> index{};
              index[(axis + 1) % 3] = other0;
              index[(axis + 2) % 3] = other1;
              std::array<const Point*, 3> line{};
              for (unsigned int position = 0; position < 3; ++position)
              {
                index[axis] = position;
                line[position] = &data.points[stencil[(index[0] * 3 + index[1]) * 3 + index[2]]];
              }
              std::array<double, 3> midpoint{};
              for (unsigned int component = 0; component < 3; ++component)
              {
                midpoint[component] = .5 * (line[0]->value[component] + line[2]->value[component]);
              }
              const auto difference = error_at(midpoint, *line[1], geometry.nominal_velocity_cm_ns);
              c.directional_error_um[axis] = std::max(c.directional_error_um[axis], std::max(difference[0], difference[1]));
            }
          }
        }
        if (c.angular_valid && std::max(c.max_transverse_error_um, c.max_time_error_um) <= policy.guard_um)
        {
          c.leaf_class = static_cast<int>(LeafClass::Certified);
          c.fallback_reason = static_cast<int>(FallbackReason::None);
          return;
        }
        c.leaf_class = static_cast<int>(LeafClass::NumericallyUnresolved);
        c.fallback_reason = static_cast<int>(c.angular_valid ? FallbackReason::Numerical : FallbackReason::AngularBranch);
        if (!c.angular_valid)
        {
          return;
        }
      }
      if (c.depth >= policy.max_depth)
      {
        c.fallback_reason = static_cast<int>(FallbackReason::DepthBudget);
        return;
      }
      const double volume = .5 * (c.hi[1] * c.hi[1] - c.lo[1] * c.lo[1]) * (c.hi[0] - c.lo[0]) * (c.hi[2] - c.lo[2]);
      pending.push({volume, id});
    };
    evaluate(root_cell);
    while (!pending.empty() && !fatal)
    {
      const auto id = pending.top().id;
      pending.pop();
      auto& c = data.cells[id];
      if (root_nodes > policy.max_nodes_per_root || policy.max_nodes_per_root - root_nodes < 2)
      {
        c.fallback_reason = static_cast<int>(FallbackReason::NodeBudget);
        continue;
      }
      const Bounds bounds{c.lo, c.hi};
      const auto axis = c.leaf_class == static_cast<int>(LeafClass::NumericallyUnresolved)
                            ? static_cast<int>(std::max_element(c.directional_error_um.begin(), c.directional_error_um.end()) - c.directional_error_um.begin())
                            : static_cast<int>(std::max_element(c.directional_reach_changes.begin(), c.directional_reach_changes.end()) - c.directional_reach_changes.begin());
      const double midpoint = .5 * (bounds.lo[axis] + bounds.hi[axis]);
      if (!(midpoint > bounds.lo[axis] && midpoint < bounds.hi[axis]))
      {
        c.fallback_reason = static_cast<int>(FallbackReason::DepthBudget);
        continue;
      }
      const auto depth = c.depth;
      if (!memory_available(data.points.size(), data.cells.size() + 2))
      {
        fatal = true;
        break;
      }
      const auto lower_id = data.cells.size();
      data.cells[id].split_axis = axis;
      data.cells[id].child_id = {{static_cast<std::int64_t>(lower_id), static_cast<std::int64_t>(lower_id + 1)}};
      // Reserve/evaluate siblings together; future splits are chosen from all
      // unresolved leaves rather than recursively spending quota in one corner.
      for (unsigned int child = 0; child < 2; ++child)
      {
        Cell cell;
        cell.id = data.cells.size();
        cell.root_id = root;
        cell.depth = depth + 1;
        const auto child_bounds = ChildBounds(bounds, axis, child);
        cell.lo = child_bounds.lo;
        cell.hi = child_bounds.hi;
        data.cells.push_back(cell);
      }
      root_nodes += 2;
      evaluate(lower_id);
      evaluate(lower_id + 1);
    }
    if (fatal)
    {
      return false;
    }
    if ((root - first) % 512 == 0)
    {
      log << "Adaptive roots completed " << root - first + 1 << '/' << end - first << " points=" << data.points.size() << " cells=" << data.cells.size() << '\n';
    }
  }
  for (const auto& cell : data.cells)
  {
    if (cell.split_axis != -1)
    {
      continue;
    }
    if (cell.leaf_class == static_cast<int>(LeafClass::Certified))
    {
      ++summary.certified_leaves;
    }
    else
    {
      ++summary.direct_leaves;
    }
    if (cell.fallback_reason >= static_cast<int>(FallbackReason::DepthBudget) && cell.fallback_reason <= static_cast<int>(FallbackReason::QueryBudget))
    {
      ++summary.budget_leaves;
    }
  }
  if (!Validate(data, error, false))
  {
    return false;
  }
  output = std::move(data);
  return true;
}
